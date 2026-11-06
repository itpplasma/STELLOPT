"""Run the native BOOTSJ edge-gradient regression when built."""
import os
import subprocess
import unittest


@unittest.skipUnless(os.environ.get("STELLOPT_BONGRID_TEST_EXE"),
                     "set STELLOPT_BONGRID_TEST_EXE to the native regression executable")
class TestBongridEdgeGradient(unittest.TestCase):
    def test_edge_gradients(self):
        result = subprocess.run([os.environ["STELLOPT_BONGRID_TEST_EXE"]],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS: BOOTSJ edge temperature gradients", result.stdout)


if __name__ == "__main__":
    unittest.main()
