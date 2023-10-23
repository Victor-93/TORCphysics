from unittest import TestCase
from TORCphysics import Environment, EnvironmentFactory, SiteFactory
from TORCphysics import binding_model as bm
from TORCphysics import unbinding_model as ubm
from TORCphysics import effect_model as em


# TODO: test the last functions
class TestEnvironment(TestCase):

    # Reads environment csv with different conditions for binding models.
    def test_environment_binding_csv(self):
        # FROM CSV
        # Cases to test csv:
        #  1. Name = None; no model.
        #  2. Name + oparams=None; Binding model with default params.
        #  3. Name + oparams; Model with params.
        #  4. Name = Model with defaults
        site_list = []
        environment_file = 'environment_binding.csv'
        csv_environment = EnvironmentFactory(filename=environment_file, site_list=site_list)
        self.assertEqual(len(csv_environment.get_environment_list()), 4)  # All loaded correctly
        self.assertEqual(csv_environment.environment_list[0].binding_model, None)  # Check specifics...
        self.assertEqual(csv_environment.environment_list[1].binding_model.k_on, 0.01)
        self.assertEqual(csv_environment.environment_list[2].binding_model.width, 0.12)
        self.assertEqual(csv_environment.environment_list[3].binding_model.width, 0.012)

    # Reads environment csv with an incorrect model name. It tests that the error is raised
    def test_environment_binding_csv_wrong_name(self):
        site_list = []

        # Check wrong model name
        environment_file2 = 'environment_binding_wrong_name.csv'
        with self.assertRaises(ValueError) as context:
            EnvironmentFactory(filename=environment_file2, site_list=site_list)
        self.assertEqual(str(context.exception), 'Could not recognise binding model Poisson')

    # Loads manually defined environmentals.
    def test_environment_binding_manual(self):
        site_list = []

        # MANUALLY DEFINED
        # Test cases for manual environmentals: #
        #  1. Model = B class.
        #  2. Model = B class + defaults
        #  3. Model = B Class + oparams
        #  4. Name + oparams dict.
        #  5. Model = No B class.

        Poisson_model = bm.PoissonBinding()

        topoI_model_default = bm.TopoIRecognition()
        oparams = {'width': 0.01, 'threshold': 0.2, 'k_on': 2.0}
        topoI_model_params = bm.TopoIRecognition(**oparams)
        NoBClass = type
        # 1. B class
        e_Poisson = Environment(e_type='RNAP', name='test1', site_list=[], concentration=0.1, size=100,
                                effective_size=50, site_type='gene', binding_model=Poisson_model)
        # 2. B class + defaults
        e_topoI_default = Environment(e_type='RNAP', name='test2', site_list=[], concentration=0.1, size=100,
                                      effective_size=50, site_type='gene', binding_model=topoI_model_default)
        # 3. B class + oparams
        e_topoI_oparams = Environment(e_type='RNAP', name='test3', site_list=[], concentration=0.1, size=100,
                                      effective_size=50, site_type='gene', binding_model=topoI_model_params)
        # 4. Name + oparams
        e_Name_oparams = Environment(e_type='RNAP', name='test4', site_list=[], concentration=0.1, size=100,
                                     effective_size=50, site_type='gene', binding_model_name='TopoIRecognition',
                                     binding_model_oparams=oparams)
        # 5. Wrong Name
        e_wrong = Environment(e_type='RNAP', name='test5', site_list=[], concentration=0.1, size=100,
                              effective_size=50, site_type='gene', binding_model=NoBClass)
        manual_environment = EnvironmentFactory(site_list=site_list)
        manual_environment.environment_list.append(e_Poisson)
        manual_environment.environment_list.append(e_topoI_default)
        manual_environment.environment_list.append(e_topoI_oparams)
        manual_environment.environment_list.append(e_Name_oparams)
        manual_environment.environment_list.append(e_wrong)
        self.assertEqual(len(manual_environment.get_environment_list()), 5)
        self.assertEqual(manual_environment.environment_list[0].binding_model_name, 'PoissonBinding')
        self.assertEqual(manual_environment.environment_list[1].binding_model.width, 0.012)
        self.assertEqual(manual_environment.environment_list[2].binding_model.width, 0.01)
        self.assertEqual(manual_environment.environment_list[3].binding_model.width, 0.01)
        self.assertEqual(manual_environment.environment_list[4].binding_model, None)

    # Reads environment csv with different conditions for effect models.
    def test_environment_effect_csv(self):
        # FROM CSV
        # Cases to test csv:
        #  1. Name = None; no model.
        #  2. Name + oparams=None; Effect model with default params.
        #  3. Name + oparams; Model with params.
        site_list = []
        environment_file = 'environment_effect.csv'
        csv_environment = EnvironmentFactory(filename=environment_file, site_list=site_list)
        self.assertEqual(len(csv_environment.get_environment_list()), 3)  # All loaded correctly
        self.assertEqual(csv_environment.environment_list[0].effect_model, None)  # Check specifics...
        self.assertEqual(csv_environment.environment_list[1].effect_model.velocity, 30)
        self.assertEqual(csv_environment.environment_list[2].effect_model.velocity, 20)

    # Reads environment csv with an incorrect model name. It tests that the error is raised. This for effect model
    def test_environment_effect_csv_wrong_name(self):
        site_list = []

        # Check wrong model name
        environment_file2 = 'environment_effect_wrong_name.csv'
        with self.assertRaises(ValueError) as context:
            EnvironmentFactory(filename=environment_file2, site_list=site_list)
        self.assertEqual(str(context.exception), 'Could not recognise effect model RNAP')

    # Loads manually defined environmentals with effect models
    def test_environment_effect_manual(self):
        site_list = []

        # MANUALLY DEFINED
        # Test cases for manual environmentals: #
        #  1. Model = E class + defaults
        #  2. Model = E Class + oparams
        #  3. Name + oparams dict.
        #  4. Model = No E class.

        RNAPUniform_default = em.RNAPUniform()
        oparams = {'velocity': 0.01, 'gamma': 0.2}
        RNAPUniform_params = em.RNAPUniform(**oparams)
        NoEClass = type
        # 1. E Model + defaults
        environment1 = Environment(e_type='RNAP', name='test1', site_list=[], concentration=0.1, size=100,
                                   effective_size=50, site_type='gene', effect_model=RNAPUniform_default)
        # 2. E class + oparams
        environment2 = Environment(e_type='RNAP', name='test2', site_list=[], concentration=0.1, size=100,
                                   effective_size=50, site_type='gene', effect_model=RNAPUniform_params)
        # 3. Name + oparams
        environment3 = Environment(e_type='RNAP', name='test3', site_list=[], concentration=0.1, size=100,
                                   effective_size=50, site_type='gene', effect_model_name='RNAPUniform',
                                   effect_model_oparams=oparams)
        # 4. Wrong class
        environment4 = Environment(e_type='RNAP', name='test5', site_list=[], concentration=0.1, size=100,
                                   effective_size=50, site_type='gene', effect_model=NoEClass)

        manual_environment = EnvironmentFactory(site_list=site_list)
        manual_environment.environment_list.append(environment1)
        manual_environment.environment_list.append(environment2)
        manual_environment.environment_list.append(environment3)
        manual_environment.environment_list.append(environment4)

        self.assertEqual(len(manual_environment.get_environment_list()), 4)  # All loaded correctly
        self.assertEqual(manual_environment.environment_list[0].effect_model.velocity, 30)
        self.assertEqual(manual_environment.environment_list[1].effect_model.velocity, 0.01)
        self.assertEqual(manual_environment.environment_list[2].effect_model.velocity, 0.01)
        self.assertEqual(manual_environment.environment_list[3].effect_model, None)

    # Reads environment csv with different conditions for unbinding models.
    def test_environment_unbinding_csv(self):
        # FROM CSV
        # Cases to test csv:
        #  1. Name = None; no model.
        #  2. Name + oparams=None; Binding model with default params.
        #  3. Name + oparams; Model with params.
        #  4. Wrong model name; should handle the situation?
        site_list = []
        environment_file = 'environment_unbinding.csv'
        csv_environment = EnvironmentFactory(filename=environment_file, site_list=site_list)
        self.assertEqual(len(csv_environment.get_environment_list()), 3)  # All loaded correctly
        self.assertEqual(csv_environment.environment_list[0].unbinding_model, None)  # Check specifics...
        self.assertEqual(csv_environment.environment_list[1].unbinding_model.k_off, 0.001)
        self.assertEqual(csv_environment.environment_list[2].unbinding_model.k_off, 2.5)

    # Reads environment csv with an incorrect model name. It tests that the error is raised
    def test_environment_unbinding_csv_wrong_name(self):
        site_list = []

        # Check wrong model name
        environment_file = 'environment_unbinding_wrong_name.csv'
        with self.assertRaises(ValueError) as context:
            EnvironmentFactory(filename=environment_file, site_list=site_list)
        self.assertEqual(str(context.exception), 'Could not recognise unbinding model Poisson')

    # Loads manually defined environmentals.
    def test_environment_unbinding_manual(self):
        site_list = []

        # MANUALLY DEFINED
        # Test cases for manual environmentals: #
        #  1. Model = UB class + defaults
        #  2. Model = UB Class + oparams
        #  3. Name + oparams dict.
        #  4. Model = No UB class.

        model_default = ubm.PoissonUnBinding()
        oparams = {'k_off': 10.0}
        model_params = ubm.PoissonUnBinding(**oparams)
        NoUBClass = type
        # 1. B class + defaults
        e_default = Environment(e_type='RNAP', name='test1', site_list=[], concentration=0.1, size=100,
                                effective_size=50, site_type='gene', unbinding_model=model_default)
        # 2. B class + oparams
        e_oparams = Environment(e_type='RNAP', name='test2', site_list=[], concentration=0.1, size=100,
                                effective_size=50, site_type='gene', unbinding_model=model_params)
        # 3. Name + oparams
        e_Name_oparams = Environment(e_type='RNAP', name='test3', site_list=[], concentration=0.1, size=100,
                                     effective_size=50, site_type='gene', unbinding_model_name='PoissonUnBinding',
                                     unbinding_model_oparams=oparams)
        # 4. Wrong Name
        e_wrong = Environment(e_type='RNAP', name='test4', site_list=[], concentration=0.1, size=100,
                              effective_size=50, site_type='gene', unbinding_model=NoUBClass)
        manual_environment = EnvironmentFactory(site_list=site_list)
        manual_environment.environment_list.append(e_default)
        manual_environment.environment_list.append(e_oparams)
        manual_environment.environment_list.append(e_Name_oparams)
        manual_environment.environment_list.append(e_wrong)
        self.assertEqual(len(manual_environment.get_environment_list()), 4)
        self.assertEqual(manual_environment.environment_list[0].unbinding_model.k_off, 0.001)
        self.assertEqual(manual_environment.environment_list[1].unbinding_model.k_off, 10.0)
        self.assertEqual(manual_environment.environment_list[2].unbinding_model.k_off, 10.0)
        self.assertEqual(manual_environment.environment_list[3].unbinding_model, None)

    # Test the incorrect case in which effective_size > size
    def test_effective_size_gt_size(self):
        with self.assertRaises(ValueError) as context:
            Environment(e_type='RNAP', name='test_effective_size', site_list=[], concentration=0.1, size=100,
                        effective_size=150, site_type='gene')
        self.assertEqual(str(context.exception), 'Error: effective_size > size')

class TestEnvironmentFactory(TestCase):

    # TODO: You need to fix these two tests because you modified your code!
    # Test the environment is not empty and that it loaded topoI correctly
    def test_EnvironmentFactory(self):
        site_list = SiteFactory("test_inputs/sites_1_gene.csv").get_site_list()
        sf = EnvironmentFactory(filename="test_inputs/environment.csv", site_list=site_list)
        self.assertGreater(len(sf.get_environment_list()), 0, "Empty environment list")
        sf_list = sf.get_environment_list()
        self.assertEqual("topoI", sf_list[0].name, "Did not load topoI correctly")

    def test_empty_environment(self):
        site_list = SiteFactory("..//sites_1_gene.csv").get_site_list()
        sf = EnvironmentFactory(filename="..//empty_environment.csv", site_list=site_list)
        self.assertEqual(len(sf.get_environment_list()), 0, "Environment not empty")
