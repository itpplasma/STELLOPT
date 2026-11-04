"""Run the native full-mesh geometry regression when built."""
import os
import subprocess
import unittest


@unittest.skipUnless(os.environ.get("STELLOPT_MAP_TEST_EXE"),
                     "set STELLOPT_MAP_TEST_EXE to the native regression executable")
class TestFullMeshMap(unittest.TestCase):
    def test_radial_intervals(self):
        result = subprocess.run([os.environ["STELLOPT_MAP_TEST_EXE"]],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS: full-mesh radial intervals", result.stdout)


if __name__ == "__main__":
    unittest.main()
