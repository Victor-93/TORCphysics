from unittest import TestCase
from TORCphysics import Enzyme, EnzymeFactory, SiteFactory


class TestEnzyme(TestCase):

    # Checks it's not empty and that it loaded the origin correctly
    def test_EnzymeFactory(self):
        site_list = SiteFactory("../sites.csv").get_site_list()
        sf = EnzymeFactory("../enzymes.csv", site_list)
        self.assertGreater(len(sf.get_enzyme_list()), 0, "Empty enzyme list")
        sf_list = sf.get_enzyme_list()
        self.assertEqual("origin", sf_list[0].enzyme_type, "Did not load origin correctly")

    def test_empty_enzyme(self):
        site_list = SiteFactory("test_inputs/sites_1_gene.csv").get_site_list()
        sf = EnzymeFactory("test_inputs/empty_enzymes.csv", site_list)
        self.assertEqual(len(sf.get_enzyme_list()), 0, "Not empty enzyme list")
