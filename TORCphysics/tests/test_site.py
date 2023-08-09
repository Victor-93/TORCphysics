from unittest import TestCase
from TORCphysics import Site, SiteFactory


class TestSite(TestCase):

    def test_SiteFactory(self):
        sf = SiteFactory("sites.csv")
        self.assertGreater(len(sf.get_site_list()), 0, "Empty site list")
        sf_list = sf.get_site_list()
        self.assertEqual("tetA", sf_list[0].name, "Did not load tetA correctly")

#    def test_set_up_oparams(self):
#        self.fail()
