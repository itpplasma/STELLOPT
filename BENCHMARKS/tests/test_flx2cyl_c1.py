import os
import subprocess
import unittest


@unittest.skipUnless(os.environ.get("STELLOPT_C1_TEST_EXE"),
                     "set STELLOPT_C1_TEST_EXE to the native regression executable")
class TestFlx2cylC1(unittest.TestCase):
    def test_c1_radial_interpolation(self):
        result = subprocess.run([os.environ["STELLOPT_C1_TEST_EXE"]],
                                capture_output=True, text=True, check=False)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS: C1 radial interpolation", result.stdout)


if __name__ == "__main__":
    unittest.main()
