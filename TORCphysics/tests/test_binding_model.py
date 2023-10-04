from unittest import TestCase
from TORCphysics import params
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em


class TestBindingModel(TestCase):

    def test_get_binding_model(self):
        # Cases:
        #  1.- b_model=None, oparams=whatever, oparams_file=whatever, model_name=None -> None, None, None, None.
        #  2.- b_model=None, oparams=dict, oparams_file=whatever, model_name = 'PoissonBinding' -> Model
        #  3.- b_model=None, oparams=None, oparams_file=None, model_name = 'PoissonBinding' -> Model
        #  4.- b_model=None, oparams=None, oparams_file=str, model_name = 'PoissonBinding' -> Model
        #  5.- b_model=BindingModel, oparams=whatever, oparams_file=whatever, model_name = whatever -> Model
        #    - The BindingModel should already be parametrised
        #  6.- b_model=EffectModel, oparams=whatever, oparams_file=whatever, model_name = whatever -> None x 4

        # Test 1
        binding_model, binding_model_name, binding_oparams_file, binding_model_oparams = (
            bm.get_binding_model(name='test1', b_model=None, model_name=None, oparams_file=None, oparams=None))
        self.assertEqual(binding_model, None)

        # Test 2
        binding_model, binding_model_name, binding_oparams_file, binding_model_oparams = (
            bm.get_binding_model(name='test2', b_model=None, model_name='PoissonBinding',
                                 oparams_file=None, oparams={'k_on': 0.2}))
        self.assertEqual(binding_model_name, 'PoissonBinding')
        self.assertEqual(binding_model.k_on, 0.2)
        self.assertEqual(binding_oparams_file, None)

        # Test 3
        binding_model, binding_model_name, binding_oparams_file, binding_model_oparams = (
            bm.get_binding_model(name='test3', b_model=None, model_name='PoissonBinding',
                                 oparams_file=None, oparams=None))
        self.assertEqual(binding_model_name, 'PoissonBinding')
        self.assertEqual(binding_model.k_on, params.k_on)

        # Test 4
        binding_model, binding_model_name, binding_oparams_file, binding_model_oparams = (
            bm.get_binding_model(name='test4', b_model=None, model_name='PoissonBinding',
                                 oparams_file='test_inputs/test_binding/PoissonBinding_params1.csv', oparams=None))
        self.assertEqual(binding_model_name, 'PoissonBinding')
        self.assertEqual(binding_model.k_on, 0.888)

        # Test 5
        my_model = bm.PoissonBinding(**{'k_on': 1.0})

        binding_model, binding_model_name, binding_oparams_file, binding_model_oparams = (
            bm.get_binding_model(name='test5', b_model=my_model, model_name='whatever',
                                 oparams_file='whatever', oparams='whatever'))
        self.assertEqual(binding_model_name, 'PoissonBinding')
        self.assertEqual(binding_model.k_on, 1.0)

        # Test 6
        effect_model = em.RNAPUniform()

        binding_model, binding_model_name, binding_oparams_file, binding_model_oparams = (
            bm.get_binding_model(name='test6', b_model=effect_model, model_name='whatever',
                                 oparams_file='whatever', oparams='whatever'))
        self.assertEqual(binding_model, None)
        self.assertEqual(binding_model_name, None)
        self.assertEqual(binding_oparams_file, None)
        self.assertEqual(binding_model_oparams, None)

    def test_assign_binding_models(self):
        # Test cases:
        #  1.- model_name='PoissonBinding' -> Model
        #  2.- model_name='PoissonBinding', filename=file -> Model
        #  3.- model_name='PoissonBinding', oparams=dict -> Model
        #  4.- model_name='PoissonBinding', filename=file, oparams=dict -> Model
        #  5.- model_name='TopoIRecognition', oparams=dict -> Model
        #  6.- model_name = 'WrongName' -> Error
        # Test 1
        my_model = bm.assign_binding_model(model_name='PoissonBinding')
        self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
        self.assertEqual(my_model.k_on, params.k_on)

        # Test 2
        my_model = bm.assign_binding_model(model_name='PoissonBinding',
                                           oparams_file='test_inputs/test_binding/PoissonBinding_params1.csv')
        self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
        self.assertEqual(my_model.k_on, 0.888)

        # Test 3
        my_model = bm.assign_binding_model(model_name='PoissonBinding',
                                           **{'k_on': 1.0})
        self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
        self.assertEqual(my_model.k_on, 1.0)

        # Test 4
        my_model = bm.assign_binding_model(model_name='PoissonBinding',
                                           oparams_file='test_inputs/test_binding/PoissonBinding_params1.csv',
                                           **{'k_on': 1.0})
        self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
        self.assertEqual(my_model.k_on, 1.0)  # DICTIONARY IS PRIORITY

        # Test 5
        my_model = bm.assign_binding_model(model_name='TopoIRecognition',
                                           **{'width': 0.01, 'threshold': 0.2, 'k_on': 2.0})
        self.assertEqual(my_model.__class__.__name__, 'TopoIRecognition')
        self.assertEqual(my_model.width, 0.01)
        self.assertEqual(my_model.threshold, 0.2)
        self.assertEqual(my_model.k_on, 2.0)

        # Test 6
        model_name = 'WrongName'
        with self.assertRaises(ValueError) as context:
            my_model = bm.assign_binding_model(model_name=model_name)
        self.assertEqual(str(context.exception), 'Could not recognise binding model ' + model_name)

# TODO: Test PoissonBinding.
