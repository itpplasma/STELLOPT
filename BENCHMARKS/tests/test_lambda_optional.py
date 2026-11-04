"""Run the native LIBSTELL optional-output regression when built."""
import os
import subprocess
import unittest


@unittest.skipUnless(os.environ.get("STELLOPT_LAMBDA_TEST_EXE"),
                     "set STELLOPT_LAMBDA_TEST_EXE to the native regression executable")
class TestLambdaOptional(unittest.TestCase):
    def test_optional_outputs(self):
        result = subprocess.run([os.environ["STELLOPT_LAMBDA_TEST_EXE"]],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("PASS: optional lambda outputs", result.stdout)


if __name__ == "__main__":
    unittest.main()
