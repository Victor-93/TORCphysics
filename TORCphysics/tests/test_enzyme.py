from unittest import TestCase
from TORCphysics import Site
from TORCphysics import Enzyme, EnzymeFactory, SiteFactory
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as ubm

# TODO: We still need to do the other effect tests and the unbinding tests
site_gene1 = Site(site_type='gene', name='test_gene1', start=100, end=500, k_on=3.00)
site_gene2 = Site(site_type='gene', name='test_gene2', start=600, end=800, k_on=3.00)
site_gene3 = Site(site_type='gene', name='test_gene3', start=1200, end=1000, k_on=3.00)
site_tf = Site(site_type='TF', name='test_TF', start=1200, end=1000, k_on=3.00)
site_list1 = [site_gene1, site_gene2, site_gene3, site_tf]


class TestEnzyme(TestCase):

    # Reads an Enzyme csv file, where the site does not exist in the site_list.
    def test_Enzyme_bad_site(self):
        enzyme_file = 'test_inputs/test_enzyme/enzyme_bad_site.csv'
        with self.assertRaises(ValueError) as context:
            EnzymeFactory(filename=enzyme_file, site_list=site_list1)
        self.assertEqual(str(context.exception), 'Error, (bound) Enzymes must be linked to a Site')

    # Tests EnzymeFactory with the possible bad inputs
    def test_bad_EnzymeFactory(self):
        enzyme_file = 'test_inputs/test_enzyme/enzyme_effect.csv'

        # Filename given but no site_list
        with self.assertRaises(ValueError) as context:
            EnzymeFactory(filename=enzyme_file)
        self.assertEqual(str(context.exception), 'Error in EnzymeFactory. filename provided but site_list is missing.')

        # site_list given but it is not a list
        with self.assertRaises(ValueError) as context:
            EnzymeFactory(site_list='Hola')
        self.assertEqual(str(context.exception), 'Error in EnzymeFactory. site_list must be a list if given.')

        # Filename given but site_list is an empty list
        with self.assertRaises(ValueError) as context:
            EnzymeFactory(filename=enzyme_file, site_list=[])
        self.assertEqual(str(context.exception), 'Error in EnzymeFactory. filename provided but empty site_list.')

    # Reads enzyme csv with different conditions for effect models.
    def test_enzyme_effect_csv(self):
        # FROM CSV
        # Cases to test csv:
        #  1. Name = None; no model.
        #  2. Name + oparams=None; Effect model with default params.
        #  3. Name + oparams; Model with params.
        enzyme_file = 'test_inputs/test_enzyme/enzyme_effect.csv'
        csv_enzyme = EnzymeFactory(filename=enzyme_file, site_list=site_list1)
        self.assertEqual(len(csv_enzyme.get_enzyme_list()), 3)  # All loaded correctly
        self.assertEqual(csv_enzyme.enzyme_list[0].effect_model, None)  # Check specifics...
        self.assertEqual(csv_enzyme.enzyme_list[1].effect_model.velocity, 30)
        self.assertEqual(csv_enzyme.enzyme_list[2].effect_model.velocity, 20)

    # Reads environment csv with an incorrect model name. It tests that the error is raised. This for effect model
    def test_enzyme_effect_csv_wrong_name(self):

        # Check wrong model name
        enzyme_file = 'test_inputs/test_enzyme/enzyme_effect_wrong_name.csv'
        with self.assertRaises(ValueError) as context:
            EnzymeFactory(filename=enzyme_file, site_list=site_list1)
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
        enzyme1 = Enzyme(e_type='RNAP', name='test1', site=site_list1[0], size=100, effective_size=50, position=30,
                         twist=0.0, superhelical=0.0, effect_model=RNAPUniform_default)
        # 2. E class + oparams
        enzyme2 = Enzyme(e_type='RNAP', name='test2', site=site_list1[1], size=100, effective_size=50, position=300,
                         twist=0.0, superhelical=0.0, effect_model=RNAPUniform_params)
        # 3. Name + oparams
        enzyme3 = Enzyme(e_type='RNAP', name='test3', site=site_list1[2], size=100, effective_size=50, position=600,
                         twist=0.0, superhelical=0.0, effect_model_name='RNAPUniform', effect_model_oparams=oparams)
        # 4. Wrong class
        enzyme4 = Enzyme(e_type='RNAP', name='test4', site=site_list1[0], size=100, effective_size=50, position=1000,
                         twist=0.0, superhelical=0.0,effect_model=NoEClass)

        manual_enzyme = EnzymeFactory()
        manual_enzyme.enzyme_list.append(enzyme1)
        manual_enzyme.enzyme_list.append(enzyme2)
        manual_enzyme.enzyme_list.append(enzyme3)
        manual_enzyme.enzyme_list.append(enzyme4)

        self.assertEqual(len(manual_enzyme.get_enzyme_list()), 4)  # All loaded correctly
        self.assertEqual(manual_enzyme.enzyme_list[0].effect_model.velocity, 30)
        self.assertEqual(manual_enzyme.enzyme_list[1].effect_model.velocity, 0.01)
        self.assertEqual(manual_enzyme.enzyme_list[2].effect_model.velocity, 0.01)
        self.assertEqual(manual_enzyme.enzyme_list[3].effect_model, None)

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
