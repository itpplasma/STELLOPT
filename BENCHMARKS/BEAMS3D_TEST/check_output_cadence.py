#!/usr/bin/env python3
"""Stored guiding-centre points must be T_END/NPOINC apart and end at T_END."""
import argparse
from pathlib import Path
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
    t_end, npoinc = 1.0e-4, 100
    with tempfile.TemporaryDirectory(prefix='beams3d-cadence-') as temporary:
        run = Path(temporary)
        text = (source/'input.ORBITS').read_text()
        # Several integration steps per stored point (dt <= 1e-7 s).
        assert text.count('40*1.0E-3') == 1 and text.count('NPOINC = 1000') == 1
        text = text.replace('40*1.0E-3', f'40*{t_end:.1E}').replace(
            'NPOINC = 1000', f'NPOINC = {npoinc}')
        (run/'input.ORBITS').write_text(text)
        shutil.copy2(source/'wout_ORBITS.nc', run/'wout_ORBITS.nc')
        result = subprocess.run([*shlex.split(args.launcher), str(executable),
            '-vmec', 'ORBITS', '-plasma'], cwd=run, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode or not (run/'beams3d_ORBITS.h5').exists():
            raise RuntimeError(result.stdout[-8000:])
        with h5py.File(run/'beams3d_ORBITS.h5') as data:
            time = data['time_lines'][:]
            radius = data['R_lines'][:]
            state = data['end_state'][:]
            last = data['t_end'][:]
        confined = np.where(state == 0)[0]
        assert len(confined) > 0, 'No confined marker to check'
        for j in confined:
            stored = time[j][radius[j] > 0]
            np.testing.assert_allclose(stored, np.linspace(0., t_end, npoinc+1),
                                       rtol=1e-9, atol=1e-15, err_msg=f'marker {j+1}')
            np.testing.assert_allclose(last[j], t_end, rtol=1e-9)
        print(f'BEAMS3D output cadence checks ({len(confined)} markers): PASS')


if __name__ == '__main__':
    main()
