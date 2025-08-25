import unittest
import os

def run_all_tests():
    # Discover and run all tests in the current directory and subdirectories
    loader = unittest.TestLoader()
    start_dir = os.path.dirname(__file__)
    suite = loader.discover(start_dir, pattern='test_*.py')

    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    # Exit with appropriate status code
    if result.wasSuccessful():
        exit(0)
    else:
        exit(1)

if __name__ == '__main__':
    run_all_tests()
