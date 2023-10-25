from unittest import TestCase
from TORCphysics import Site, Enzyme, Environment, params, utils
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
import numpy as np

# Let's just define a site_list to use for various tests
site_gene1 = Site(site_type='gene', name='test_gene1', start=100, end=500, k_on=3.00)
site_gene2 = Site(site_type='gene', name='test_gene2', start=600, end=800, k_on=3.00)
site_gene3 = Site(site_type='gene', name='test_gene3', start=1200, end=1000, k_on=3.00)
site_tf = Site(site_type='TF', name='test_TF', start=1200, end=1000, k_on=3.00)
site_list1 = [site_gene1, site_gene2, site_gene3, site_tf]

# And list of enzymes (RNAPs without effect model)
enzyme1 = Enzyme(e_type='RNAP', name='test1', site=site_list1[0], size=30, effective_size=15, position=30,
                 twist=0.0, superhelical=0.0)
enzyme2 = Enzyme(e_type='RNAP', name='test2', site=site_list1[1], size=30, effective_size=15, position=300,
                 twist=0.0, superhelical=0.0)
enzyme3 = Enzyme(e_type='RNAP', name='test3', site=site_list1[2], size=30, effective_size=15, position=600,
                 twist=0.0, superhelical=0.0)
enzyme4 = Enzyme(e_type='RNAP', name='test4', site=site_list1[0], size=30, effective_size=15, position=1000,
                 twist=0.0, superhelical=0.0)
enzyme_list1 = [enzyme1, enzyme2, enzyme3, enzyme4]


class TestEffectModel(TestCase):
    def test_get_effect_model(self):
        # Cases:
        #  1.- e_model=None, oparams=whatever, oparams_file=whatever, model_name=None -> None, None, None, None.
        #  2.- e_model=None, oparams=dict, oparams_file=whatever, model_name = 'RNAPUniform' -> Model
        #  3.- e_model=None, oparams=None, oparams_file=None, model_name = 'RNAPUniform' -> Model
        #  4.- e_model=None, oparams=None, oparams_file=str, model_name = 'RNAPUniform' -> Model
        #  5.- e_model=BindingModel, oparams=whatever, oparams_file=whatever, model_name = whatever -> Model
        #    - The BindingModel should already be parametrised
        #  6.- e_model=EffectModel, oparams=whatever, oparams_file=whatever, model_name = whatever -> None x 4

        # Test 1
        effect_model, effect_model_name, effect_oparams_file, effect_model_oparams = (
            em.get_effect_model(name='test1', e_model=None, model_name=None, oparams_file=None, oparams=None))
        self.assertEqual(effect_model, None)

        # Test 2
        effect_model, effect_model_name, effect_oparams_file, effect_model_oparams = (
            em.get_effect_model(name='test2', e_model=None, model_name='RNAPUniform',
                                oparams_file=None, oparams={'velocity': 30.0, 'gamma': 0.02}))
        self.assertEqual(effect_model_name, 'RNAPUniform')
        self.assertEqual(effect_model.velocity, 30.0)
        self.assertEqual(effect_model.gamma, 0.02)
        self.assertEqual(effect_oparams_file, None)

        # Test 3
        effect_model, effect_model_name, effect_oparams_file, effect_model_oparams = (
            em.get_effect_model(name='test3', e_model=None, model_name='RNAPUniform',
                                oparams_file=None, oparams=None))
        self.assertEqual(effect_model_name, 'RNAPUniform')
        self.assertEqual(effect_model.velocity, params.v0)
        self.assertEqual(effect_model.gamma, params.gamma)

        # Test 4
        effect_model, effect_model_name, effect_oparams_file, effect_model_oparams = (
            em.get_effect_model(name='test4', e_model=None, model_name='RNAPUniform',
                                oparams_file='RNAPUniform_params.csv', oparams=None))
        self.assertEqual(effect_model_name, 'RNAPUniform')
        self.assertEqual(effect_model.velocity, 20.0)
        self.assertEqual(effect_model.gamma, 2.0)

        # Test 5
        my_model = em.RNAPUniform(**{'velocity': 30.0, 'gamma': 0.02})

        effect_model, effect_model_name, effect_oparams_file, effect_model_oparams = (
            em.get_effect_model(name='test5', e_model=my_model, model_name='whatever',
                                oparams_file='whatever', oparams='whatever'))
        self.assertEqual(effect_model_name, 'RNAPUniform')
        self.assertEqual(effect_model.velocity, 30.0)
        self.assertEqual(effect_model.gamma, 0.02)

        # Test 6
        binding_model = bm.PoissonBinding()

        effect_model, effect_model_name, effect_oparams_file, effect_model_oparams = (
            em.get_effect_model(name='test6', e_model=binding_model, model_name='whatever',
                                oparams_file='whatever', oparams='whatever'))
        self.assertEqual(effect_model, None)
        self.assertEqual(effect_model_name, None)
        self.assertEqual(effect_oparams_file, None)
        self.assertEqual(effect_model_oparams, None)

    def test_assign_effect_models(self):
        # Test cases:
        #  1.- model_name='RNAPUniform' -> Model
        #  2.- model_name='RNAPUniform', filename=file -> Model
        #  3.- model_name='RNAPUniform', oparams=dict -> Model
        #  4.- model_name='RNAPUniform', filename=file, oparams=dict -> Model
        #  5.- model_name = 'WrongName' -> Error

        # Test 1
        my_model = em.assign_effect_model(model_name='RNAPUniform')
        self.assertEqual(my_model.__class__.__name__, 'RNAPUniform')
        self.assertEqual(my_model.velocity, params.v0)
        self.assertEqual(my_model.gamma, params.gamma)

        # Test 2
        my_model = em.assign_effect_model(model_name='RNAPUniform',
                                          oparams_file='RNAPUniform_params.csv')
        self.assertEqual(my_model.__class__.__name__, 'RNAPUniform')
        self.assertEqual(my_model.velocity, 20.0)
        self.assertEqual(my_model.gamma, 2.0)

        # Test 3
        my_model = em.assign_effect_model(model_name='RNAPUniform',
                                          **{'velocity': 30.0, 'gamma': 0.02})
        self.assertEqual(my_model.__class__.__name__, 'RNAPUniform')
        self.assertEqual(my_model.velocity, 30.0)
        self.assertEqual(my_model.gamma, 0.02)

        # Test 4
        my_model = em.assign_effect_model(model_name='RNAPUniform',
                                          oparams_file='RNAPUniform_params.csv',
                                          **{'velocity': 30.0, 'gamma': 0.02})
        self.assertEqual(my_model.__class__.__name__, 'RNAPUniform')
        self.assertEqual(my_model.velocity, 30.0)  # DICTIONARY IS PRIORITY
        self.assertEqual(my_model.gamma, 0.02)

        # Test 5
        model_name = 'WrongName'
        with self.assertRaises(ValueError) as context:
            my_model = em.assign_effect_model(model_name=model_name)
        self.assertEqual(str(context.exception), 'Could not recognise effect model ' + model_name)


class TestRNAPUniform(TestCase):

    def test_RNAPUniform(self):
        # For each test case, we should have the RNAPUniform Model with params, and get an Effect.
        #  The test cases should work for a various timesteps.
        #  Test cases:
        #  1.- filename=None, no oparams -> Model & Effect
        #  2.- filename=file, no oparams
        #  3.- filename=None, oparams
        #  4.- filename=file, oparams

        for dt in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:
            # Test 1
            # ---------------------------
            my_model = em.RNAPUniform()  # default
            my_enzyme = Enzyme(e_type='RNAP', name='RNAPUniform_test', site=site_list1[0], size=30,
                               effective_size=15, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.
            # Note that the complete code should do this on its own, but since this is a test,
            # I'm actually doing it manually

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = my_enzyme.direction * my_model.velocity * dt
            correct_twist = my_model.gamma * my_model.velocity * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'RNAPUniform')
            self.assertEqual(my_enzyme.effect_model.velocity, params.v0)
            self.assertEqual(my_enzyme.effect_model.gamma, params.gamma)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 2
            # ---------------------------
            my_model = em.RNAPUniform(filename='RNAPUniform_params.csv')
            my_enzyme = Enzyme(e_type='RNAP', name='RNAPUniform_test', site=site_list1[0], size=30,
                               effective_size=15, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = my_enzyme.direction * my_model.velocity * dt
            correct_twist = my_model.gamma * my_model.velocity * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'RNAPUniform')
            self.assertEqual(my_model.velocity, 20.0)
            self.assertEqual(my_model.gamma, 2.0)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 3
            # ---------------------------
            my_model = em.RNAPUniform(**{'velocity': 30.0, 'gamma': 0.02})
            my_enzyme = Enzyme(e_type='RNAP', name='RNAPUniform_test', site=site_list1[0], size=30,
                               effective_size=15, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = my_enzyme.direction * my_model.velocity * dt
            correct_twist = my_model.gamma * my_model.velocity * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'RNAPUniform')
            self.assertEqual(my_model.velocity, 30.0)
            self.assertEqual(my_model.gamma, 0.02)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 4
            # ---------------------------
            my_model = em.RNAPUniform(filename='RNAPUniform_params.csv',
                                      **{'velocity': 30.0, 'gamma': 0.02})  # oparams is priority!
            my_enzyme = Enzyme(e_type='RNAP', name='RNAPUniform_test', site=site_list1[0], size=30,
                               effective_size=15, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = my_enzyme.direction * my_model.velocity * dt
            correct_twist = my_model.gamma * my_model.velocity * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'RNAPUniform')
            self.assertEqual(my_model.velocity, 30.0)
            self.assertEqual(my_model.gamma, 0.02)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

    def test_RNAPUniform_bad_csv(self):
        # test 1 = bad velocity
        filename = 'RNAPUniform_badvelocity.csv'
        with self.assertRaises(ValueError) as context:
            em.RNAPUniform(filename=filename)
        self.assertEqual(str(context.exception), 'Error, velocity parameter missing in csv file for RNAPUniform')

        # test 2 = bad gamma
        filename = 'RNAPUniform_badgamma.csv'
        with self.assertRaises(ValueError) as context:
            em.RNAPUniform(filename=filename)
        self.assertEqual(str(context.exception), 'Error, gamma parameter missing in csv file for RNAPUniform')


class TestTopoIUniform(TestCase):

    def test_TopoIUniform(self):
        # For each test case, we should have the TopoIUniform Model with params, and get an Effect.
        #  The test cases should work for a various timesteps.
        #  Test cases:
        #  1.- filename=None, no oparams -> Model & Effect
        #  2.- filename=file, no oparams
        #  3.- filename=None, oparams
        #  4.- filename=file, oparams

        for dt in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:
            # Test 1
            # ---------------------------
            my_model = em.TopoIUniform()  # default
            my_enzyme = Enzyme(e_type='RNAP', name='TopoIUniform_test', site=site_list1[0], size=100,
                               effective_size=80, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.
            # Note that the complete code should do this on its own, but since this is a test,
            # I'm actually doing it manually

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = 0.0
            correct_twist = 0.5 * my_model.k_cat * params.w0 * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'TopoIUniform')
            self.assertEqual(my_enzyme.effect_model.k_cat, params.topoI_uniform_k_cat)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 2
            # ---------------------------
            my_model = em.TopoIUniform(filename='TopoIUniform_params.csv')
            my_enzyme = Enzyme(e_type='RNAP', name='TopoIUniform_test', site=site_list1[0], size=100,
                               effective_size=80, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.

            # Calculate true values
            correct_position = 0.0
            correct_twist = 0.5 * my_model.k_cat * params.w0 * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'TopoIUniform')
            self.assertEqual(my_enzyme.effect_model.k_cat, 10.2)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 3
            # ---------------------------
            my_model = em.TopoIUniform(**{'k_cat': 13.0})
            my_enzyme = Enzyme(e_type='RNAP', name='TopoIUniform_test', site=site_list1[0], size=100,
                               effective_size=80, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = 0.0
            correct_twist = 0.5 * my_model.k_cat * params.w0 * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'TopoIUniform')
            self.assertEqual(my_enzyme.effect_model.k_cat, 13.0)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 4
            # ---------------------------
            my_model = em.TopoIUniform(filename='TopoIUniform_params.csv',
                                       **{'k_cat': 13.0})  # oparams is priority!
            my_enzyme = Enzyme(e_type='RNAP', name='TopoIUniform_test', site=site_list1[0], size=100,
                               effective_size=80, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = 0.0
            correct_twist = 0.5 * my_model.k_cat * params.w0 * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'TopoIUniform')
            self.assertEqual(my_enzyme.effect_model.k_cat, 13.0)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

    def test_TopoIUniform_bad_csv(self):
        # test 1 = bad velocity
        filename = 'RNAPUniform_badvelocity.csv'
        with self.assertRaises(ValueError) as context:
            em.TopoIUniform(filename=filename)
        self.assertEqual(str(context.exception), 'Error, k_cat parameter missing in csv file for TopoIUniform')


class TestGyraseUniform(TestCase):

    def test_GyraseUniform(self):
        # For each test case, we should have the GyraseUniform Model with params, and get an Effect.
        #  The test cases should work for a various timesteps.
        #  Test cases:
        #  1.- filename=None, no oparams -> Model & Effect
        #  2.- filename=file, no oparams
        #  3.- filename=None, oparams
        #  4.- filename=file, oparams

        for dt in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:
            # Test 1
            # ---------------------------
            my_model = em.GyraseUniform()  # default
            my_enzyme = Enzyme(e_type='RNAP', name='GyraseUniform_test', site=site_list1[0], size=100,
                               effective_size=80, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.
            # Note that the complete code should do this on its own, but since this is a test,
            # I'm actually doing it manually

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = 0.0
            correct_twist = 0.5 * my_model.k_cat * params.w0 * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'GyraseUniform')
            self.assertEqual(my_enzyme.effect_model.k_cat, params.gyra_uniform_k_cat)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - abs(correct_twist), 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - abs(correct_twist), 0.000001)

            # Test 2
            # ---------------------------
            my_model = em.GyraseUniform(filename='GyraseUniform_params.csv')
            my_enzyme = Enzyme(e_type='RNAP', name='GyraseUniform_test', site=site_list1[0], size=100,
                               effective_size=80, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = 0.0
            correct_twist = 0.5 * my_model.k_cat * params.w0 * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'GyraseUniform')
            self.assertEqual(my_enzyme.effect_model.k_cat, -10.2)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - abs(correct_twist), 0.0001)
            self.assertLessEqual(abs(my_effect.twist_right) - abs(correct_twist), 0.0001)

            # Test 3
            # ---------------------------
            my_model = em.GyraseUniform(**{'k_cat': -13.0})
            my_enzyme = Enzyme(e_type='RNAP', name='GyraseUniform_test', site=site_list1[0], size=100,
                               effective_size=80, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = 0.0
            correct_twist = 0.5 * my_model.k_cat * params.w0 * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'GyraseUniform')
            self.assertEqual(my_enzyme.effect_model.k_cat, -13.0)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - abs(correct_twist), 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - abs(correct_twist), 0.000001)

            # Test 4
            # ---------------------------
            my_model = em.GyraseUniform(filename='GyraseUniform_params.csv',
                                        **{'k_cat': -13.0})  # oparams is priority!
            my_enzyme = Enzyme(e_type='RNAP', name='GyraseUniform_test', site=site_list1[0], size=100,
                               effective_size=80, position=500, twist=0.0, superhelical=0.0, effect_model=my_model)
            my_list = [enzyme1, enzyme2, my_enzyme, enzyme3, enzyme4]
            my_index = 2  # According to my enzyme position, the index is 2.

            # Calculate effect
            my_effect = my_enzyme.effect_model.calculate_effect(my_index, my_enzyme, my_list, dt)

            # Calculate true values
            correct_position = 0.0
            correct_twist = 0.5 * my_model.k_cat * params.w0 * dt

            # And test
            self.assertEqual(my_enzyme.effect_model.__class__.__name__, 'GyraseUniform')
            self.assertEqual(my_enzyme.effect_model.k_cat, -13.0)

            self.assertEqual(my_effect.position, correct_position)
            self.assertLessEqual(abs(my_effect.twist_left) - abs(correct_twist), 0.000001)
            self.assertLessEqual(abs(my_effect.twist_right) - abs(correct_twist), 0.000001)

    def test_GyraseUniform_bad_csv(self):
        # test 1 = bad velocity
        filename = 'RNAPUniform_badvelocity.csv'
        with self.assertRaises(ValueError) as context:
            em.GyraseUniform(filename=filename)
        self.assertEqual(str(context.exception), 'Error, k_cat parameter missing in csv file for GyraseUniform')


class TestTopoIContinuum(TestCase):

    def test_TopoIContinuum(self):
        # For each test case, we should have the TopoIContinuum Model with params, and get an Effect.
        #  The test cases should work for a various timesteps.
        #  Test cases:
        #  1.- filename=None, no oparams -> Model & Effect
        #  2.- filename=file, no oparams
        #  3.- filename=None, oparams
        #  4.- filename=file, oparams

        for dt in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:
            # Test 1
            # ---------------------------
            my_model = em.TopoIContinuum()  # default
            my_list = [enzyme1, enzyme2, enzyme3, enzyme4]
            # Note that the complete code should do this on its own, but since this is a test,
            # I'm actually doing it manually
            my_environment = Environment(e_type='topo', name='topoIC_test', site_list=[], concentration=0.1, size=0.0,
                                         effective_size=0.0, site_type='none', effect_model=my_model)
            for i, enzyme in enumerate(my_list):
                if i == 3:  # Let's not touch the end.
                    continue
                # Calculate effect and add it to the list
                my_effect = my_environment.effect_model.calculate_effect(concentration=my_environment.concentration,
                                                                         index=i, z=enzyme, z_list=my_list, dt=dt)
                # Calculate true values
                enzyme_n = [e for e in my_list if e.position > enzyme.position][0]  # Enzyme on the right
                a = my_environment.concentration * my_model.k_cat * dt
                try:
                    b = 1 + np.exp((enzyme.superhelical - my_model.threshold) / my_model.width)
                    supercoiling_removed = a / b
                except OverflowError as oe:
                    supercoiling_removed = 0.0

                correct_position = 0.0
                correct_twist = utils.calculate_twist_from_sigma(enzyme, enzyme_n, supercoiling_removed)

                # And test
                self.assertEqual(my_environment.effect_model.__class__.__name__, 'TopoIContinuum')
                self.assertEqual(my_environment.effect_model.k_cat, params.topo_sam_kcat)
                self.assertEqual(my_environment.effect_model.threshold, params.topo_sam_threshold)
                self.assertEqual(my_environment.effect_model.width, params.topo_sam_width)

                self.assertEqual(my_effect.position, correct_position)
                self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
                self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 2
            # ---------------------------
            my_model = em.TopoIContinuum(filename='TopoIContinuum_params.csv')
            my_environment = Environment(e_type='topo', name='topoIC_test', site_list=[], concentration=0.1, size=0.0,
                                         effective_size=0.0, site_type='none', effect_model=my_model)
            for i, enzyme in enumerate(my_list):
                if i == 3:  # Let's not touch the end.
                    continue
                # Calculate effect and add it to the list
                my_effect = my_environment.effect_model.calculate_effect(concentration=my_environment.concentration,
                                                                         index=i, z=enzyme, z_list=my_list, dt=dt)
                # Calculate true values
                enzyme_n = [e for e in my_list if e.position > enzyme.position][0]  # Enzyme on the right
                a = my_environment.concentration * my_model.k_cat * dt
                try:
                    b = 1 + np.exp((enzyme.superhelical - my_model.threshold) / my_model.width)
                    supercoiling_removed = a / b
                except OverflowError as oe:
                    supercoiling_removed = 0.0

                correct_position = 0.0
                correct_twist = utils.calculate_twist_from_sigma(enzyme, enzyme_n, supercoiling_removed)

                # And test
                self.assertEqual(my_environment.effect_model.__class__.__name__, 'TopoIContinuum')
                self.assertEqual(my_environment.effect_model.k_cat, 0.003)
                self.assertEqual(my_environment.effect_model.threshold, -0.06)
                self.assertEqual(my_environment.effect_model.width, 0.02)

                self.assertEqual(my_effect.position, correct_position)
                self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
                self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 3
            # ---------------------------
            my_model = em.TopoIContinuum(**{'k_cat': 0.005, 'width': 0.03, 'threshold': -0.01})
            my_environment = Environment(e_type='topo', name='topoIC_test', site_list=[], concentration=0.1, size=0.0,
                                         effective_size=0.0, site_type='none', effect_model=my_model)
            for i, enzyme in enumerate(my_list):
                if i == 3:  # Let's not touch the end.
                    continue
                # Calculate effect and add it to the list
                my_effect = my_environment.effect_model.calculate_effect(concentration=my_environment.concentration,
                                                                         index=i, z=enzyme, z_list=my_list, dt=dt)
                # Calculate true values
                enzyme_n = [e for e in my_list if e.position > enzyme.position][0]  # Enzyme on the right
                a = my_environment.concentration * my_model.k_cat * dt
                try:
                    b = 1 + np.exp((enzyme.superhelical - my_model.threshold) / my_model.width)
                    supercoiling_removed = a / b
                except OverflowError as oe:
                    supercoiling_removed = 0.0

                correct_position = 0.0
                correct_twist = utils.calculate_twist_from_sigma(enzyme, enzyme_n, supercoiling_removed)

                # And test
                self.assertEqual(my_environment.effect_model.__class__.__name__, 'TopoIContinuum')
                self.assertEqual(my_environment.effect_model.k_cat, 0.005)
                self.assertEqual(my_environment.effect_model.threshold, -0.01)
                self.assertEqual(my_environment.effect_model.width, 0.03)

                self.assertEqual(my_effect.position, correct_position)
                self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
                self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 4
            # ---------------------------
            my_model = em.TopoIContinuum(filename='TopoIContinuum_params.csv',
                                         **{'k_cat': 0.005, 'width': 0.03, 'threshold': -0.01})  # oparams is priority!
            my_environment = Environment(e_type='topo', name='topoIC_test', site_list=[], concentration=0.1, size=0.0,
                                         effective_size=0.0, site_type='none', effect_model=my_model)
            for i, enzyme in enumerate(my_list):
                if i == 3:  # Let's not touch the end.
                    continue
                # Calculate effect and add it to the list
                my_effect = my_environment.effect_model.calculate_effect(concentration=my_environment.concentration,
                                                                         index=i, z=enzyme, z_list=my_list, dt=dt)
                # Calculate true values
                enzyme_n = [e for e in my_list if e.position > enzyme.position][0]  # Enzyme on the right
                a = my_environment.concentration * my_model.k_cat * dt
                try:
                    b = 1 + np.exp((enzyme.superhelical - my_model.threshold) / my_model.width)
                    supercoiling_removed = a / b
                except OverflowError as oe:
                    supercoiling_removed = 0.0

                correct_position = 0.0
                correct_twist = utils.calculate_twist_from_sigma(enzyme, enzyme_n, supercoiling_removed)

                # And test
                self.assertEqual(my_environment.effect_model.__class__.__name__, 'TopoIContinuum')
                self.assertEqual(my_environment.effect_model.k_cat, 0.005)
                self.assertEqual(my_environment.effect_model.threshold, -0.01)
                self.assertEqual(my_environment.effect_model.width, 0.03)

                self.assertEqual(my_effect.position, correct_position)
                self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
                self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

    def test_TopoIContinuum_bad_csv(self):
        # test 1 = bad k_cat
        filename = 'TopoIContinuum_badkcat.csv'
        with self.assertRaises(ValueError) as context:
            em.TopoIContinuum(filename=filename)
        self.assertEqual(str(context.exception), 'Error, k_cat parameter missing in csv file for TopoIContinuum')

        # test 2 = bad threshold
        filename = 'TopoIContinuum_badthreshold.csv'
        with self.assertRaises(ValueError) as context:
            em.TopoIContinuum(filename=filename)
        self.assertEqual(str(context.exception), 'Error, threshold parameter missing in csv file for TopoIContinuum')

        # test 3 = bad width
        filename = 'TopoIContinuum_badwidth.csv'
        with self.assertRaises(ValueError) as context:
            em.TopoIContinuum(filename=filename)
        self.assertEqual(str(context.exception), 'Error, width parameter missing in csv file for TopoIContinuum')


class TestGyraseContinuum(TestCase):

    def test_GyraseContinuum(self):
        # For each test case, we should have the GyraseContinuum Model with params, and get an Effect.
        #  The test cases should work for a various timesteps.
        #  Test cases:
        #  1.- filename=None, no oparams -> Model & Effect
        #  2.- filename=file, no oparams
        #  3.- filename=None, oparams
        #  4.- filename=file, oparams

        for dt in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:
            # Test 1
            # ---------------------------
            my_model = em.GyraseContinuum()  # default
            my_list = [enzyme1, enzyme2, enzyme3, enzyme4]
            # Note that the complete code should do this on its own, but since this is a test,
            # I'm actually doing it manually
            my_environment = Environment(e_type='topo', name='GyraseC_test', site_list=[], concentration=0.1, size=0.0,
                                         effective_size=0.0, site_type='none', effect_model=my_model)
            for i, enzyme in enumerate(my_list):
                if i == 3:  # Let's not touch the end.
                    continue
                # Calculate effect and add it to the list
                my_effect = my_environment.effect_model.calculate_effect(concentration=my_environment.concentration,
                                                                         index=i, z=enzyme, z_list=my_list, dt=dt)
                # Calculate true values
                enzyme_n = [e for e in my_list if e.position > enzyme.position][0]  # Enzyme on the right
                a = my_environment.concentration * my_model.k_cat * dt
                try:
                    b = 1 + np.exp((enzyme.superhelical - my_model.threshold) / my_model.width)
                    supercoiling_removed = a / b
                except OverflowError as oe:
                    supercoiling_removed = 0.0

                correct_position = 0.0
                correct_twist = utils.calculate_twist_from_sigma(enzyme, enzyme_n, supercoiling_removed)

                # And test
                self.assertEqual(my_environment.effect_model.__class__.__name__, 'GyraseContinuum')
                self.assertEqual(my_environment.effect_model.k_cat, params.gyra_sam_kcat)
                self.assertEqual(my_environment.effect_model.threshold, params.gyra_sam_threshold)
                self.assertEqual(my_environment.effect_model.width, params.gyra_sam_width)

                self.assertEqual(my_effect.position, correct_position)
                self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
                self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 2
            # ---------------------------
            my_model = em.GyraseContinuum(filename='GyraseContinuum_params.csv')
            my_environment = Environment(e_type='topo', name='GyraseC_test', site_list=[], concentration=0.1, size=0.0,
                                         effective_size=0.0, site_type='none', effect_model=my_model)
            for i, enzyme in enumerate(my_list):
                if i == 3:  # Let's not touch the end.
                    continue
                # Calculate effect and add it to the list
                my_effect = my_environment.effect_model.calculate_effect(concentration=my_environment.concentration,
                                                                         index=i, z=enzyme, z_list=my_list, dt=dt)
                # Calculate true values
                enzyme_n = [e for e in my_list if e.position > enzyme.position][0]  # Enzyme on the right
                a = my_environment.concentration * my_model.k_cat * dt
                try:
                    b = 1 + np.exp((enzyme.superhelical - my_model.threshold) / my_model.width)
                    supercoiling_removed = a / b
                except OverflowError as oe:
                    supercoiling_removed = 0.0

                correct_position = 0.0
                correct_twist = utils.calculate_twist_from_sigma(enzyme, enzyme_n, supercoiling_removed)

                # And test
                self.assertEqual(my_environment.effect_model.__class__.__name__, 'GyraseContinuum')
                self.assertEqual(my_environment.effect_model.k_cat, 0.003)
                self.assertEqual(my_environment.effect_model.threshold, 0.06)
                self.assertEqual(my_environment.effect_model.width, 0.02)

                self.assertEqual(my_effect.position, correct_position)
                self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
                self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 3
            # ---------------------------
            my_model = em.GyraseContinuum(**{'k_cat': 0.005, 'width': 0.03, 'threshold': 0.01})
            my_environment = Environment(e_type='topo', name='GyraseC_test', site_list=[], concentration=0.1, size=0.0,
                                         effective_size=0.0, site_type='none', effect_model=my_model)
            for i, enzyme in enumerate(my_list):
                if i == 3:  # Let's not touch the end.
                    continue
                # Calculate effect and add it to the list
                my_effect = my_environment.effect_model.calculate_effect(concentration=my_environment.concentration,
                                                                         index=i, z=enzyme, z_list=my_list, dt=dt)
                # Calculate true values
                enzyme_n = [e for e in my_list if e.position > enzyme.position][0]  # Enzyme on the right
                a = my_environment.concentration * my_model.k_cat * dt
                try:
                    b = 1 + np.exp((enzyme.superhelical - my_model.threshold) / my_model.width)
                    supercoiling_removed = a / b
                except OverflowError as oe:
                    supercoiling_removed = 0.0

                correct_position = 0.0
                correct_twist = utils.calculate_twist_from_sigma(enzyme, enzyme_n, supercoiling_removed)

                # And test
                self.assertEqual(my_environment.effect_model.__class__.__name__, 'GyraseContinuum')
                self.assertEqual(my_environment.effect_model.k_cat, 0.005)
                self.assertEqual(my_environment.effect_model.threshold, 0.01)
                self.assertEqual(my_environment.effect_model.width, 0.03)

                self.assertEqual(my_effect.position, correct_position)
                self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
                self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

            # Test 4
            # ---------------------------
            my_model = em.GyraseContinuum(filename='GyraseContinuum_params.csv',
                                         **{'k_cat': 0.005, 'width': 0.03, 'threshold': 0.01})  # oparams is priority!
            my_environment = Environment(e_type='topo', name='GyraseC_test', site_list=[], concentration=0.1, size=0.0,
                                         effective_size=0.0, site_type='none', effect_model=my_model)
            for i, enzyme in enumerate(my_list):
                if i == 3:  # Let's not touch the end.
                    continue
                # Calculate effect and add it to the list
                my_effect = my_environment.effect_model.calculate_effect(concentration=my_environment.concentration,
                                                                         index=i, z=enzyme, z_list=my_list, dt=dt)
                # Calculate true values
                enzyme_n = [e for e in my_list if e.position > enzyme.position][0]  # Enzyme on the right
                a = my_environment.concentration * my_model.k_cat * dt
                try:
                    b = 1 + np.exp((enzyme.superhelical - my_model.threshold) / my_model.width)
                    supercoiling_removed = a / b
                except OverflowError as oe:
                    supercoiling_removed = 0.0

                correct_position = 0.0
                correct_twist = utils.calculate_twist_from_sigma(enzyme, enzyme_n, supercoiling_removed)

                # And test
                self.assertEqual(my_environment.effect_model.__class__.__name__, 'GyraseContinuum')
                self.assertEqual(my_environment.effect_model.k_cat, 0.005)
                self.assertEqual(my_environment.effect_model.threshold, 0.01)
                self.assertEqual(my_environment.effect_model.width, 0.03)

                self.assertEqual(my_effect.position, correct_position)
                self.assertLessEqual(abs(my_effect.twist_left) - correct_twist, 0.000001)
                self.assertLessEqual(abs(my_effect.twist_right) - correct_twist, 0.000001)

    def test_GyraseContinuum_bad_csv(self):
        # test 1 = bad k_cat
        filename = 'TopoIContinuum_badkcat.csv'
        with self.assertRaises(ValueError) as context:
            em.GyraseContinuum(filename=filename)
        self.assertEqual(str(context.exception), 'Error, k_cat parameter missing in csv file for GyraseContinuum')

        # test 2 = bad threshold
        filename = 'TopoIContinuum_badthreshold.csv'
        with self.assertRaises(ValueError) as context:
            em.GyraseContinuum(filename=filename)
        self.assertEqual(str(context.exception), 'Error, threshold parameter missing in csv file for GyraseContinuum')

        # test 3 = bad width
        filename = 'TopoIContinuum_badwidth.csv'
        with self.assertRaises(ValueError) as context:
            em.GyraseContinuum(filename=filename)
        self.assertEqual(str(context.exception), 'Error, width parameter missing in csv file for GyraseContinuum')


class TestEffect(TestCase):

    # Just test that it assigns the data correctly.
    def test_effect_1(self):
        for index, enzyme in enumerate(enzyme_list1):
            my_effect = em.Effect(index=index, position=enzyme.position,
                                  twist_left=enzyme.twist, twist_right=enzyme.twist)
            self.assertEqual(my_effect.index, index)
            self.assertEqual(my_effect.position, enzyme.position)
            self.assertEqual(my_effect.twist_left, enzyme.twist)
            self.assertEqual(my_effect.twist_right, enzyme.twist)


class TestEffectFunctions(TestCase):

    def test_get_start_end_c(self):
        nbp = 5
        my_site = Site(s_type='test_site', name='u1', start=1, end=nbp, k_min=0, k_max=0, s_model_name=None,
                       oparams=None)
        my_enzyme = Enzyme(e_type='test', name='YO', site=my_site, position=2, size=1, k_cat=0.0, k_off=0.0,
                           twist=0, superhelical=0)
        position_left, position_right = utils.get_start_end_c(my_enzyme, my_enzyme, nbp)
        self.assertEqual(-2, position_left, "Incorrectly calculated the left position")
        self.assertEqual(7, position_right, "Incorrectly calculated the right position")

    def test_calculate_length(self):
        nbp = 2000
        my_site = Site(s_type='test_site', name='u1', start=1, end=nbp, k_min=0, k_max=0, s_model_name=None,
                       oparams=None)
        enzyme1 = Enzyme(e_type='test', name='YO', site=my_site, position=200, size=50, k_cat=0.0, k_off=0.0,
                         twist=0, superhelical=0)
        enzyme2 = Enzyme(e_type='test', name='YO', site=my_site, position=1200, size=50, k_cat=0.0, k_off=0.0,
                         twist=0, superhelical=0)
        length = utils.calculate_length(enzyme1, enzyme2)
        self.assertEqual(length, 950, "Incorrectly calculated the length")
