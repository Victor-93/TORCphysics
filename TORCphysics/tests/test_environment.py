from unittest import TestCase
from TORCphysics import Environment, EnvironmentFactory, SiteFactory


class TestEnvironment(TestCase):

    # Test the environment is not empty and that it loaded topoI correctly
    def test_EnvironmentFactory(self):
        site_list = SiteFactory("test_inputs/sites_1_gene.csv").get_site_list()
        sf = EnvironmentFactory("test_inputs/environment.csv", site_list)
        self.assertGreater(len(sf.get_environment_list()), 0, "Empty environment list")
        sf_list = sf.get_environment_list()
        self.assertEqual("topoI", sf_list[0].name, "Did not load topoI correctly")

    def test_empty_environment(self):
        site_list = SiteFactory("test_inputs/sites_1_gene.csv").get_site_list()
        sf = EnvironmentFactory("test_inputs/empty_environment.csv", site_list)
        self.assertEqual(len(sf.get_environment_list()), 0, "Environment not empty")
