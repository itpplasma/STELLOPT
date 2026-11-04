#!/usr/bin/env python3
"""The default toroidal grid must close at two pi in double precision."""
import argparse
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--launcher', default='')
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    executable = args.executable.resolve()
    with tempfile.TemporaryDirectory(prefix='beams3d-period-') as temporary:
        run = Path(temporary)
        text = (source/'input.ORBITS_loss').read_text()
        text, count = re.subn(r'(?im)^\s*PHIMAX\s*=.*$', '', text)
        assert count == 1, 'Expected one PHIMAX in the published orbit fixture'
        text = text.replace('40*1.0E-3', '40*1.0E-9')
        (run/'input.ORBITS_loss').write_text(text)
        shutil.copy2(source/'wout_ORBITS_loss.nc', run/'wout_ORBITS_loss.nc')
        result = subprocess.run([*shlex.split(args.launcher), str(executable),
            '-vmec', 'ORBITS_loss', '-plasma'], cwd=run, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode or not (run/'beams3d_ORBITS_loss.h5').exists():
            raise RuntimeError(result.stdout[-8000:])
        with h5py.File(run/'beams3d_ORBITS_loss.h5') as data:
            phi = data['phiaxis'][:]
        np.testing.assert_allclose(phi, np.linspace(0.,2*np.pi,len(phi)),
                                   rtol=0., atol=2e-14)
        print('BEAMS3D default period checks: PASS')


if __name__ == '__main__':
    main()
