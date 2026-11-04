"""Run the native analytic-flux regression when built."""
import os
import subprocess
import unittest


@unittest.skipUnless(os.environ.get("STELLOPT_FLUX_TEST_EXE"),
                     "set STELLOPT_FLUX_TEST_EXE to the native regression executable")
class TestVectorPotentialFlux(unittest.TestCase):
    def test_flux_normalization(self):
        result = subprocess.run([os.environ["STELLOPT_FLUX_TEST_EXE"]],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS: vector potential flux normalization", result.stdout)


if __name__ == "__main__":
    unittest.main()
