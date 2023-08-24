from unittest import TestCase
from TORCphysics import Environment, EnvironmentFactory, SiteFactory
from TORCphysics import binding_model as bm


class TestEnvironment(TestCase):

    def test_model_given(self):

        model = bm.PoissonBinding()
        environmental = Environment(e_type='RNAP', name='test', site_list=[], concentration=0.1, size=100, eff_size=50,
                                    site_type='gene', binding_model=model)
#        self.assertGreater(len(environmental.get_environment_list()), 0, "Empty environment list")

        #e_type, name, site_list, concentration, size, eff_size, site_type,
#                 binding_model_name=None, binding_oparams_file=None,
#                 effect_model_name=None, effect_oparams_file=None,
#                 unbinding_model_name=None, unbinding_oparams_file=None,
#                 binding_model=None, effect_model=None, unbinding_model=None):

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
