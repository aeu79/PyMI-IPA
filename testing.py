import unittest


"""
    Unit testing for the MI_IPA python implementation.

    Usage:
    import testing
    test(all)    #Performs all the tests
    testx()     # Test x
    profiler()   # Makes a profile with the sample data.

    TODO:
"""


def profiler():
    import cProfile
    import mi_ipa_main
    cProfile.run('mi_ipa_main.main()', 'profile_stats','tottime')

class TestMiipa(unittest.TestCase):

    def test_1(self):
        self.assertEqual(sum([1, 2, 3]), 6, "Should be 6")


    def test_2(self):
        self.assertEqual(sum((1, 2, 2)), 5, "Should be 6")



# if __name__ == "__main__":
#     testall()

if __name__ == '__main__':
    unittest.main()