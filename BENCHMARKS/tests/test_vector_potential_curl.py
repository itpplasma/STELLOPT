"""Run the native curl(A) = B regression when built."""
import os
import subprocess
import unittest
from pathlib import Path

WOUT = Path(__file__).resolve().parents[1] / "BEAMS3D_TEST" / "wout_ORBITS.nc"


@unittest.skipUnless(os.environ.get("STELLOPT_CURL_TEST_EXE"),
                     "set STELLOPT_CURL_TEST_EXE to the native regression executable")
class TestVectorPotentialCurl(unittest.TestCase):
    def test_curl_matches_field(self):
        result = subprocess.run([os.environ["STELLOPT_CURL_TEST_EXE"], str(WOUT)],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS: vector potential curl", result.stdout)


if __name__ == "__main__":
    unittest.main()
