from unittest import TestCase
from TORCphysics import params, Environment
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em

environmental1 = Environment(e_type='RNAP', name='test1', site_list=[], concentration=0.1, size=100,
                             effective_size=50, site_type='gene')


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
                                 oparams_file='PoissonBinding_params1.csv', oparams=None))
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
                                           oparams_file='PoissonBinding_params1.csv')
        self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
        self.assertEqual(my_model.k_on, 0.888)

        # Test 3
        my_model = bm.assign_binding_model(model_name='PoissonBinding',
                                           **{'k_on': 1.0})
        self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
        self.assertEqual(my_model.k_on, 1.0)

        # Test 4
        my_model = bm.assign_binding_model(model_name='PoissonBinding',
                                           oparams_file='PoissonBinding_params1.csv',
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


class TestPoissonBinding(TestCase):

    def test_PoissonBinding(self):
        # For each test case, we should have the PoissonBinding Model with params, and 0 <= probability <=1.
        #  The test cases should work for a various timesteps.
        #  Test cases:
        #  1.- filename=None, no oparams -> Model & 0 <= probability <=1
        #  2.- filename=file, no oparams
        #  3.- filename=None, oparams
        #  4.- filename=file, oparams

        for dt in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:  #
            supercoiling = -0.06
            # Test 1
            my_model = bm.PoissonBinding()
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
            self.assertEqual(my_model.k_on, params.k_on)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 2
            my_model = bm.PoissonBinding(filename='PoissonBinding_params1.csv')
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
            self.assertEqual(my_model.k_on, 0.888)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 3
            my_model = bm.PoissonBinding(**{'k_on': 1.0})
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
            self.assertEqual(my_model.k_on, 1.0)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 4
            my_model = bm.PoissonBinding(filename='PoissonBinding_params1.csv',
                                         **{'k_on': 1.0})  # oparams is priority!
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'PoissonBinding')
            self.assertEqual(my_model.k_on, 1.0)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

    def test_PoissonBinding_bad_csv(self):
        # bad k_on
        filename = 'PoissonBinding_bad.csv'
        with self.assertRaises(ValueError) as context:
            bm.PoissonBinding(filename=filename)
        self.assertEqual(str(context.exception), 'Error, k_on parameter missing in csv file for PoissonBinding')


class TestTopoIRecognition(TestCase):

    def test_TopoIRecognition(self):
        # For each test case, we should have the TopoIRecognition Model with params, and 0 <= probability <=1.
        #  The test cases should work for a various timesteps.
        #  Test cases:
        #  1.- filename=None, no oparams -> Model & 0 <= probability <=1
        #  2.- filename=file, no oparams
        #  3.- filename=None, oparams
        #  4.- filename=file, oparams

        for dt in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:  #
            supercoiling = -0.06
            # Test 1
            my_model = bm.TopoIRecognition()
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'TopoIRecognition')
            self.assertEqual(my_model.k_on, params.topo_b_k_on)
            self.assertEqual(my_model.width, params.topo_b_w)
            self.assertEqual(my_model.threshold, params.topo_b_t)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 2
            my_model = bm.TopoIRecognition(filename='TopoIRecognition_params1.csv')
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'TopoIRecognition')
            self.assertEqual(my_model.k_on, 0.005)
            self.assertEqual(my_model.width, 0.03)
            self.assertEqual(my_model.threshold, 0.01)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 3
            my_model = bm.TopoIRecognition(**{'k_on': 1.0, 'width': 0.3, 'threshold': 0.2})

            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'TopoIRecognition')
            self.assertEqual(my_model.k_on, 1.0)
            self.assertEqual(my_model.width, 0.3)
            self.assertEqual(my_model.threshold, 0.2)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 4
            my_model = bm.TopoIRecognition(filename='TopoIRecognition_params1.csv',
                                           **{'k_on': 1.0, 'width': 0.3, 'threshold': 0.2})  # oparams is priority!
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'TopoIRecognition')
            self.assertEqual(my_model.k_on, 1.0)
            self.assertEqual(my_model.width, 0.3)
            self.assertEqual(my_model.threshold, 0.2)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

    def test_TopoIRecognition_bad_csv(self):
        # bad k_on
        filename = 'TopoIRecognition_bad.csv'
        with self.assertRaises(ValueError) as context:
            bm.TopoIRecognition(filename=filename)
        self.assertEqual(str(context.exception), 'Error, k_on parameter missing in csv file for TopoIRecognition')


class TestGyraseRecognition(TestCase):

    def test_GyraseRecognition(self):
        # For each test case, we should have the GyraseRecognition Model with params, and 0 <= probability <=1.
        #  The test cases should work for a various timesteps.
        #  Test cases:
        #  1.- filename=None, no oparams -> Model & 0 <= probability <=1
        #  2.- filename=file, no oparams
        #  3.- filename=None, oparams
        #  4.- filename=file, oparams

        for dt in [0.01, 0.1, 0.5, 1.0, 1.5, 2.0]:  #
            supercoiling = -0.06
            # Test 1
            my_model = bm.GyraseRecognition()
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'GyraseRecognition')
            self.assertEqual(my_model.k_on, params.gyra_b_k_on)
            self.assertEqual(my_model.width, params.gyra_b_w)
            self.assertEqual(my_model.threshold, params.gyra_b_t)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 2
            my_model = bm.GyraseRecognition(filename='GyraseRecognition_params1.csv')
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'GyraseRecognition')
            self.assertEqual(my_model.k_on, 0.005)
            self.assertEqual(my_model.width, 0.03)
            self.assertEqual(my_model.threshold, 0.01)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 3
            my_model = bm.GyraseRecognition(**{'k_on': 1.0, 'width': 0.3, 'threshold': 0.2})

            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'GyraseRecognition')
            self.assertEqual(my_model.k_on, 1.0)
            self.assertEqual(my_model.width, 0.3)
            self.assertEqual(my_model.threshold, 0.2)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

            # Test 4
            my_model = bm.GyraseRecognition(filename='GyraseRecognition_params1.csv',
                                            **{'k_on': 1.0, 'width': 0.3, 'threshold': 0.2})  # oparams is priority!
            probability = my_model.binding_probability(environmental=environmental1, superhelical=supercoiling, dt=dt)
            self.assertEqual(my_model.__class__.__name__, 'GyraseRecognition')
            self.assertEqual(my_model.k_on, 1.0)
            self.assertEqual(my_model.width, 0.3)
            self.assertEqual(my_model.threshold, 0.2)
            self.assertGreaterEqual(probability, 0.0)
            self.assertLessEqual(probability, 1.0)

    def test_GyraseRecognition_bad_csv(self):
        # bad k_on
        filename = 'GyraseRecognition_bad.csv'
        with self.assertRaises(ValueError) as context:
            bm.GyraseRecognition(filename=filename)
        self.assertEqual(str(context.exception), 'Error, k_on parameter missing in csv file for GyraseRecognition')
