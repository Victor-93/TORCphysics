from unittest import TestCase
from TORCphysics import params
from TORCphysics import unbinding_model as ubm
from TORCphysics import effect_model as em


class TestBindingModel(TestCase):

    # TODO: Testea esto
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
                                    oparams_file='test_inputs/test_unbinding/PoissonUnBinding_params1.csv',
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
                                              oparams_file='test_inputs/test_unbinding/PoissonUnBinding_params1.csv')
        self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
        self.assertEqual(my_model.k_off, 2.5)

        # Test 3
        my_model = ubm.assign_unbinding_model(model_name='PoissonUnBinding',
                                              **{'k_off': 1.0})
        self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
        self.assertEqual(my_model.k_off, 1.0)

        # Test 4
        my_model = ubm.assign_unbinding_model(model_name='PoissonUnBinding',
                                              oparams_file='test_inputs/test_unbinding/PoissonUnBinding_params1.csv',
                                              **{'k_off': 1.0})
        self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
        self.assertEqual(my_model.k_off, 1.0)  # DICTIONARY IS PRIORITY

        # Test 5
        model_name = 'WrongName'
        with self.assertRaises(ValueError) as context:
            my_model = ubm.assign_unbinding_model(model_name=model_name)
        self.assertEqual(str(context.exception), 'Could not recognise unbinding model ' + model_name)

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
            probability = my_model.unbinding_probability(dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
            self.assertEqual(my_model.k_off, params.k_off)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 2
            my_model = ubm.PoissonUnBinding(filename='test_inputs/test_unbinding/PoissonUnBinding_params1.csv')
            probability = my_model.unbinding_probability(dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
            self.assertEqual(my_model.k_off, 2.5)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 3
            my_model = ubm.PoissonUnBinding(**{'k_off': 1.0})
            probability = my_model.unbinding_probability(dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
            self.assertEqual(my_model.k_off, 1.0)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 4
            my_model = ubm.PoissonUnBinding(filename='test_inputs/test_unbinding/PoissonUnBinding_params1.csv',
                                            **{'k_off': 1.0})  # oparams is priority!
            probability = my_model.unbinding_probability(dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonUnBinding')
            self.assertEqual(my_model.k_off, 1.0)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)
