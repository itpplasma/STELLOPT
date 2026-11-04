"""Run the native missing-chi WOUT regression when built."""
import os
import subprocess
import unittest
from pathlib import Path

FIXTURES = Path(__file__).resolve().parents[1] / "BEAMS3D_TEST"


@unittest.skipUnless(os.environ.get("STELLOPT_CHI_TEST_EXE"),
                     "set STELLOPT_CHI_TEST_EXE to the native regression executable")
class TestWoutMissingChi(unittest.TestCase):
    def test_reconstructed_chi(self):
        result = subprocess.run([os.environ["STELLOPT_CHI_TEST_EXE"],
                                 str(FIXTURES / "wout_ORBITS.nc"),
                                 str(FIXTURES / "wout_ORBITS_nochi.nc")],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS: missing chi reconstruction", result.stdout)


if __name__ == "__main__":
    unittest.main()
