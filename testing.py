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

    def test_csv_species_list(self):
        """Comparison of imported species from .mat and csv files."""
        # from matlab data to df
        import scipy.io as spio
        import pandas as pd
        # Read the matlab file:
        self.mat = spio.loadmat('testing/SpeciesNumbering_Standard_HKRR_dataset.mat', squeeze_me=True)
        self.col1 = list(range(len(self.mat['SpeciesNumbering_extr'])))
        self.col2 = [x[1] for x in self.mat['SpeciesNumbering_extr']]
        self.SpeciesNumbering_mat = pd.DataFrame([*zip(self.col1, self.col2)])
        self.SpeciesNumbering_mat = self.SpeciesNumbering_mat.values.tolist() # Panda's DataFrame to list
        # Now import the csv file:
        self.species_list = 'testing/SpeciesNumbering_Standard_HKRR_dataset.csv'
        self.SpeciesNumbering_csv = pd.read_csv(self.species_list, header=None, delimiter=',')
        self.SpeciesNumbering_csv = self.SpeciesNumbering_csv.values.tolist()
        #Only_unique_species_1col = self.SpeciesNumbering_csv[1].unique().tolist()
        self.assertListEqual(self.SpeciesNumbering_mat, self.SpeciesNumbering_csv, "The list of species should be the same from .mat or .csv.")

    def test_2(self):
        self.assertEqual(sum((1, 2, 2)), 5, "Should be 6")


if __name__ == '__main__':
    unittest.main(verbosity=2)

