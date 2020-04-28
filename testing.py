"""
    Unit testing for the MI_IPA python implementation.

    Usage:
    python3 testing.py      # Performs all the tests

    Individual tests:
    import testing
    testing.test_x()      # Test x
    
    Profiler:
    import testing
    testing.profiler()   # Makes a profile with the sample data (default).
"""

import unittest


# noinspection SpellCheckingInspection
def profiler(output_stats='profile_stats'):
    """
    Usage:
        import testing
        import mi_ipa_main
        testing.profiler(arguments, output_stats) # arguments not implemented yet as cProfiler crashes with variables.
    """
    import cProfile
    # import mi_ipa_main
    # cProfile.run(mi_ipa_main.main(), 'profile_stats', 'tottime') # Using this import and the following command crashes cProfile (report?):
    import pstats

    cProfile.run('mi_ipa_main.main()', 'profile_stats', 'tottime')
    profile_results = pstats.Stats(output_stats).strip_dirs()
    profile_results.sort_stats('tottime').print_stats(10)
    print("You can now obtain a graph using: gprof2dot -f pstats profile_stats | dot -Tsvg -o " + output_stats + ".svg")


class TestMiIpa(unittest.TestCase):

    def test_1(self):
        self.assertEqual(sum([1, 2, 3]), 6, "Should be 6")

    def test_2(self):
        self.assertEqual(sum((1, 2, 2)), 5, "Should be 6")


if __name__ == '__main__':
    unittest.main(verbosity=2)
