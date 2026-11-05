import os
import subprocess
import unittest


@unittest.skipUnless(os.environ.get("STELLOPT_COULOMB_TEST_EXE"),
                     "set STELLOPT_COULOMB_TEST_EXE to the native regression executable")
class TestCoulombRates(unittest.TestCase):
    def test_rates_match_nrl(self):
        result = subprocess.run([os.environ["STELLOPT_COULOMB_TEST_EXE"]],
                                capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS: BEAMS3D Coulomb rates", result.stdout)


if __name__ == "__main__":
    unittest.main()
