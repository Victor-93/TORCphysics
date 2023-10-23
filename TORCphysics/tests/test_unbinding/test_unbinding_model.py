from unittest import TestCase
from TORCphysics import params, Site, Enzyme
from TORCphysics import unbinding_model as ubm
from TORCphysics import effect_model as em

site_gene1 = Site(site_type='gene', name='test_gene1', start=100, end=500, k_on=3.00)
site_gene2 = Site(site_type='gene', name='test_gene2', start=600, end=800, k_on=3.00)
site_gene3 = Site(site_type='gene', name='test_gene3', start=1200, end=1000, k_on=3.00)
site_tf = Site(site_type='TF', name='test_TF', start=1200, end=1000, k_on=3.00)
site_list1 = [site_gene1, site_gene2, site_gene3, site_tf]

# 1. E Model + defaults
enzyme1 = Enzyme(e_type='RNAP', name='test1', site=site_list1[0], size=100, effective_size=50, position=30,
                 twist=0.0, superhelical=0.0)


class TestBindingModel(TestCase):

    def test_get_unbinding_model(self):
        # Cases:
        #  1.- ub_model=None, oparams=whatever, oparams_file=whatever, model_name=None -> None, None, None, None.
        #  2.- ub_model=None, oparams=dict, oparams_file=whatever, model_name = 'PoissonUnBinding' -> Model
        #  3.- ub_model=None, oparams=None, oparams_file=None, model_name = 'PoissonUnBinding' -> Model
        #  4.- ub_model=None, oparams=None, oparams_file=str, model_name = 'PoissonUnBinding' -> Model
        #  5.- ub_model=BindingModel, oparams=whatever, oparams_file=whatever, model_name = whatever -> Model
        #    - The BindingModel should already be parametrised
        #  6.- ub_model=EffectModel, oparams=whatever, oparams_file=whatever, model_name = whatever -> None x 4

        # Test 1
        unbinding_model, unbinding_model_name, unbinding_oparams_file, unbinding_model_oparams = (
            ubm.get_unbinding_model(name='test1', ub_model=None, model_name=None, oparams_file=None, oparams=None))
        self.assertEqual(unbinding_model, None)

        # Test 2
        unbinding_model, unbinding_model_name, unbinding_oparams_file, unbinding_model_oparams = (
            ubm.get_unbinding_model(name='test2', ub_model=None, model_name='PoissonUnBinding',
                                    oparams_file=None, oparams={'k_off': 0.2}))
        self.assertEqual(unbinding_model_name, 'PoissonUnBinding')
        self.assertEqual(unbinding_model.k_off, 0.2)
        self.assertEqual(unbinding_oparams_file, None)

        # Test 3
        unbinding_model, unbinding_model_name, unbinding_oparams_file, unbinding_model_oparams = (
            ubm.get_unbinding_model(name='test3', ub_model=None, model_name='PoissonUnBinding',
                                    oparams_file=None, oparams=None))
        self.assertEqual(unbinding_model_name, 'PoissonUnBinding')
        self.assertEqual(unbinding_model.k_off, params.k_off)

        # Test 4
        unbinding_model, unbinding_model_name, unbinding_oparams_file, unbinding_model_oparams = (
            ubm.get_unbinding_model(name='test4', ub_model=None, model_name='PoissonUnBinding',
                                    oparams_file='PoissonUnBinding_params1.csv',
                                    oparams=None))
        self.assertEqual(unbinding_model_name, 'PoissonUnBinding')
        self.assertEqual(unbinding_model.k_off, 2.5)

        # Test 5
        my_model = ubm.PoissonUnBinding(**{'k_off': 1.0})

        unbinding_model, unbinding_model_name, unbinding_oparams_file, unbinding_model_oparams = (
            ubm.get_unbinding_model(name='test5', ub_model=my_model, model_name='whatever',
                                    oparams_file='whatever', oparams='whatever'))
        self.assertEqual(unbinding_model_name, 'PoissonUnBinding')
        self.assertEqual(unbinding_model.k_off, 1.0)

        # Test 6
        effect_model = em.RNAPUniform()

        unbinding_model, unbinding_model_name, unbinding_oparams_file, unbinding_model_oparams = (
            ubm.get_unbinding_model(name='test6', ub_model=effect_model, model_name='whatever',
                                    oparams_file='whatever', oparams='whatever'))
        self.assertEqual(unbinding_model, None)
        self.assertEqual(unbinding_model_name, None)
        self.assertEqual(unbinding_oparams_file, None)
        self.assertEqual(unbinding_model_oparams, None)

    def test_assign_unbinding_models(self):
        # Test cases:
        #  1.- model_name='PoissonUnBinding' -> Model
        #  2.- model_name='PoissonUnBinding', filename=file -> Model
        #  3.- model_name='PoissonUnBinding', oparams=dict -> Model
        #  4.- model_name='PoissonUnBinding', filename=file, oparams=dict -> Model
        #  5.- model_name = 'WrongName' -> Error
        # Test 1
        my_model = ubm.assign_unbinding_model(model_name='PoissonUnBinding')
        self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
        self.assertEqual(my_model.k_off, params.k_off)

        # Test 2
        my_model = ubm.assign_unbinding_model(model_name='PoissonUnBinding',
                                              oparams_file='PoissonUnBinding_params1.csv')
        self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
        self.assertEqual(my_model.k_off, 2.5)

        # Test 3
        my_model = ubm.assign_unbinding_model(model_name='PoissonUnBinding',
                                              **{'k_off': 1.0})
        self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
        self.assertEqual(my_model.k_off, 1.0)

        # Test 4
        my_model = ubm.assign_unbinding_model(model_name='PoissonUnBinding',
                                              oparams_file='PoissonUnBinding_params1.csv',
                                              **{'k_off': 1.0})
        self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
        self.assertEqual(my_model.k_off, 1.0)  # DICTIONARY IS PRIORITY

        # Test 5
        model_name = 'WrongName'
        with self.assertRaises(ValueError) as context:
            my_model = ubm.assign_unbinding_model(model_name=model_name)
        self.assertEqual(str(context.exception), 'Could not recognise unbinding model ' + model_name)


class TestPoissonUnBindingModel(TestCase):

    def test_PoissonUnBinding(self):
        # For each test case, we should have the PoissonUnBinding Model with params, and 0 <= probability <=1.
        #  The test cases should work for a various timesteps.
        #  Test cases:
        #  1.- filename=None, no oparams -> Model & 0 <= probability <=1
        #  2.- filename=file, no oparams
        #  3.- filename=None, oparams
        #  4.- filename=file, oparams

        for dt in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:
            # Test 1
            my_model = ubm.PoissonUnBinding()
            probability = my_model.unbinding_probability(enzyme1, dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
            self.assertEqual(my_model.k_off, params.k_off)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 2
            my_model = ubm.PoissonUnBinding(filename='PoissonUnBinding_params1.csv')
            probability = my_model.unbinding_probability(enzyme1, dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
            self.assertEqual(my_model.k_off, 2.5)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 3
            my_model = ubm.PoissonUnBinding(**{'k_off': 1.0})
            probability = my_model.unbinding_probability(enzyme1, dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
            self.assertEqual(my_model.k_off, 1.0)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 4
            my_model = ubm.PoissonUnBinding(filename='PoissonUnBinding_params1.csv',
                                            **{'k_off': 1.0})  # oparams is priority!
            probability = my_model.unbinding_probability(enzyme1, dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
            self.assertEqual(my_model.k_off, 1.0)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

    def test_PoissonUnBinding_bad_csv(self):
        # Poisson Binding model that has k_on but not k_off
        filename = '..//test_binding/PoissonBinding_params1.csv'
        with self.assertRaises(ValueError) as context:
            ubm.PoissonUnBinding(filename=filename)
        self.assertEqual(str(context.exception), 'Error, k_off parameter missing in csv file for PoissonUnBinding')


class TestRNAPSimpleUnbindingModel(TestCase):

    def test_RNAPSimpleUnbinding_load(self):
        # Here we test if the model is loaded correctly. This model does not need any parameters and if they're given
        # those should be ignored.
        #  Test cases:
        #  1.- filename=None, no oparams, and checks name
        #  2.- filename=file, no oparams, and checks name
        #  3.- filename=None, oparams and checks name
        #  4.- filename=file, oparams
        my_model = ubm.RNAPSimpleUnbinding()
        self.assertEqual(my_model.__class__.__name__, 'RNAPSimpleUnbinding')

        # Test 2
        my_model = ubm.RNAPSimpleUnbinding(filename='RNAPSimpleUnbinding_params1.csv')  # This file doesn't even exist
        self.assertEqual(my_model.__class__.__name__, 'RNAPSimpleUnbinding')

        # Test 3
        my_model = ubm.RNAPSimpleUnbinding(**{'k_off': 1.0})
        self.assertEqual(my_model.__class__.__name__, 'RNAPSimpleUnbinding')

        # Test 4
        my_model = ubm.RNAPSimpleUnbinding(filename='RNAPSimpleUnbinding_params1.csv',
                                           **{'k_off': 1.0})  # oparams is priority!
        self.assertEqual(my_model.__class__.__name__, 'RNAPSimpleUnbinding')

    def test_RNAPSimpleUnbinding_probabilities(self):

        # For each test case, we should have the RNAPSimpleUnbinding Model with params, and 0 <= probability <=1.
        #  Test cases:
        #   R = RNAP, S = Start, E = End
        #   1.- S____R____E; direction +1, p =0
        #   2.- S________RE; p=1
        #   3.- S_________E____R; p=1
        #   4.- R____S_____E;  p=0
        #   5.- S____R______E; but direction=0, error?
        #   6.- E____R____S; direction -1, p =0
        #   7.- ER________S; p=1
        #   8.- R___E_________S; p=1
        #   9.- E_____S______R;  p=0
        #   10.- E____R______S; but direction=0, error?

        unbinding_model = ubm.RNAPSimpleUnbinding()
        dt = 1.0

        # This site will be used for tests 1-5
        # --------------------------------------------------------------
        test_site = Site(site_type='gene', name='test_gene', start=200, end=1200, k_on=3.00)

        # Test 1 : S____R____E; direction +1, p =0
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=500,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        probability = test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(test_enzyme.unbinding_model.__class__.__name__, 'RNAPSimpleUnbinding')
        self.assertGreaterEqual(probability, 0.0)
        self.assertLessEqual(probability, 1.0)
        self.assertEqual(probability, 0.0)

        # Test 2 : S________RE; p=1
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=1200,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        probability = test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(test_enzyme.unbinding_model.__class__.__name__, 'RNAPSimpleUnbinding')
        self.assertGreaterEqual(probability, 0.0)
        self.assertLessEqual(probability, 1.0)
        self.assertEqual(probability, 1.0)

        # Test 3 : S_________E____R; p=1
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=1500,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        probability = test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(test_enzyme.unbinding_model.__class__.__name__, 'RNAPSimpleUnbinding')
        self.assertGreaterEqual(probability, 0.0)
        self.assertLessEqual(probability, 1.0)
        self.assertEqual(probability, 1.0)

        # Test 4 : R____S_____E;  p=0
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=10,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        probability = test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(test_enzyme.unbinding_model.__class__.__name__, 'RNAPSimpleUnbinding')
        self.assertGreaterEqual(probability, 0.0)
        self.assertLessEqual(probability, 1.0)
        self.assertEqual(probability, 0.0)

        # Test 5 : S____R______E; but direction=0, error?
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=500,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        test_enzyme.direction = 0
        with self.assertRaises(ValueError) as context:
            test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(str(context.exception), 'Error. Enzyme with invalid direction in RNAPSimpleUnbinding.')

        # This site will be used for tests 6-10
        # --------------------------------------------------------------
        test_site = Site(site_type='gene', name='test_gene', start=1200, end=200, k_on=3.00)

        # Test 6 : E____R____S; direction -1, p =0
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=500,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        probability = test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(test_enzyme.unbinding_model.__class__.__name__, 'RNAPSimpleUnbinding')
        self.assertGreaterEqual(probability, 0.0)
        self.assertLessEqual(probability, 1.0)
        self.assertEqual(probability, 0.0)

        # Test 7 : ER________S; p=1
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=200,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        probability = test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(test_enzyme.unbinding_model.__class__.__name__, 'RNAPSimpleUnbinding')
        self.assertGreaterEqual(probability, 0.0)
        self.assertLessEqual(probability, 1.0)
        self.assertEqual(probability, 1.0)

        # Test 8 : R___E_________S; p=1
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=20,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        probability = test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(test_enzyme.unbinding_model.__class__.__name__, 'RNAPSimpleUnbinding')
        self.assertGreaterEqual(probability, 0.0)
        self.assertLessEqual(probability, 1.0)
        self.assertEqual(probability, 1.0)
        # Test 9 : E_____S______R;  p=0
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=1500,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        probability = test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(test_enzyme.unbinding_model.__class__.__name__, 'RNAPSimpleUnbinding')
        self.assertGreaterEqual(probability, 0.0)
        self.assertLessEqual(probability, 1.0)
        self.assertEqual(probability, 0.0)

        # Test 10 : E____R______S; but direction=0, error?
        test_enzyme = Enzyme(e_type='RNAP', name='test', site=test_site, size=50, effective_size=30, position=500,
                             twist=0.0, superhelical=0.0, unbinding_model=unbinding_model)
        test_enzyme.direction = 0
        with self.assertRaises(ValueError) as context:
            test_enzyme.unbinding_model.unbinding_probability(test_enzyme, dt)
        self.assertEqual(str(context.exception), 'Error. Enzyme with invalid direction in RNAPSimpleUnbinding.')
