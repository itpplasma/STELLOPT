#!/usr/bin/env python3

import argparse
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

import h5py
import numpy as np


def run_case(source, executable, launcher, seed):
    with tempfile.TemporaryDirectory(prefix=f"beams3d-seed-{seed}-") as temporary:
        run_dir = Path(temporary)
        text = (source / "input.ORBITS_slow").read_text()
        anchor = "  NPOINC = 1000\n"
        if text.count(anchor) != 1:
            raise ValueError("unexpected NPOINC assignment")
        text = text.replace(anchor, anchor + f"  RANDOM_SEED_IN = {seed}\n")
        (run_dir / "input.ORBITS_slow").write_text(text)
        shutil.copy2(source / "wout_ORBITS_slow.nc", run_dir / "wout_ORBITS_slow.nc")
        command = [*launcher, str(executable), "-vmec", "ORBITS_slow", "-plasma",
                   "-collisions"]
        result = subprocess.run(command, cwd=run_dir, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            output = "\n".join((result.stdout + result.stderr).splitlines()[-80:])
            raise RuntimeError(f"seed {seed} failed with status {result.returncode}\n{output}")
        data = {}
        with h5py.File(run_dir / "beams3d_ORBITS_slow.h5", "r") as handle:
            handle.visititems(lambda name, item: data.__setitem__(name, item[()])
                              if isinstance(item, h5py.Dataset) else None)
        return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--launcher", default="")
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    launcher = shlex.split(args.launcher)
    executable = args.executable.resolve()
    first = run_case(source, executable, launcher, 12345)
    second = run_case(source, executable, launcher, 12345)
    other = run_case(source, executable, launcher, 54321)
    assert first.keys() == second.keys()
    differ = [name for name in first if name != "dist_prof"
              and not np.array_equal(first[name], second[name])]
    # Ranks add to the shared distribution without a lock, so the order of
    # the sums, and with it the last bits, varies between runs.
    scale = np.max(np.abs(first["dist_prof"]))
    if not np.allclose(first["dist_prof"], second["dist_prof"],
                       rtol=0.0, atol=1.0e-12*scale):
        differ.append("dist_prof")
    if differ:
        raise AssertionError(f"same seed, different datasets: {differ}")
    for name in ("end_state", "epower_prof", "ipower_prof"):
        assert name in first, name
    if np.array_equal(first["R_lines"], other["R_lines"]):
        raise AssertionError("different seeds gave identical orbits")
    print("BEAMS3D random seed checks: PASS")


if __name__ == "__main__":
    main()
