from unittest import TestCase
from TORCphysics import Site, SiteFactory
from TORCphysics import binding_model as bm


class TestSite(TestCase):

    # Reads site csv with different conditions for binding models.
    def test_Site_csv(self):
        # FROM CSV
        # Cases to test csv:
        #  1. No model.
        #  2. Name ; Binding model with default params.
        #  3. Name + oparams + k_on=0.0; Model with params and k_on = 0.0
        site_file = 'test_inputs/test_site/site1.csv'
        csv_site = SiteFactory(filename=site_file)
        self.assertEqual(len(csv_site.site_list), 3)  # All loaded correctly
        self.assertEqual(csv_site.site_list[0].binding_model, None)  # Check specifics...
        self.assertEqual(csv_site.site_list[1].binding_model.k_on, 100.0)
        self.assertEqual(csv_site.site_list[1].binding_model_name, 'PoissonBinding')
        self.assertEqual(csv_site.site_list[2].binding_model_name, 'TopoIRecognition')
        self.assertEqual(csv_site.site_list[2].binding_model.k_on, 0.0)

    # Reads site csv with an incorrect model name. It tests that the error is raised
    def test_environment_binding_csv_wrong_name(self):
        # Check wrong model name
        site_file = 'test_inputs/test_site/site_model_wrong_name.csv'
        with self.assertRaises(ValueError) as context:
            SiteFactory(filename=site_file)
        self.assertEqual(str(context.exception), 'Could not recognise binding model Poisson')

    # Loads manually defined sites.
    def test_sites_binding_manual(self):
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
        NoBClass = bm.PoissonUnBinding
        # 1. B class
        s_Poisson = Site(site_type='gene', name='Poisson_default',
                         start=10, end=20, k_on=3.0, binding_model=Poisson_model)
        # 2. B class + defaults
        s_topoI_default = Site(site_type='gene', name='topoI_default',
                               start=10, end=20, k_on=3.0, binding_model=topoI_model_default)
        # 3. B class + oparams
        s_topoI_oparams = Site(site_type='gene', name='topoI_oparams',
                               start=10, end=20, k_on=3.0, binding_model=topoI_model_params)
        # 4. Name + oparams
        s_name_oparams = Site(site_type='gene', name='name_oparams',
                              start=10, end=20, k_on=3.0,
                              binding_model_name='TopoIRecognition', binding_model_oparams=oparams)
        # 5. Wrong Name
        s_wrong = Site(site_type='gene', name='wrong_B_class',
                       start=10, end=20, k_on=3.0, binding_model=NoBClass)
        manual_site = SiteFactory()
        manual_site.site_list.append(s_Poisson)
        manual_site.site_list.append(s_topoI_default)
        manual_site.site_list.append(s_topoI_oparams)
        manual_site.site_list.append(s_name_oparams)
        manual_site.site_list.append(s_wrong)
        self.assertEqual(len(manual_site.site_list), 5)
        self.assertEqual(manual_site.site_list[0].binding_model_name, 'PoissonBinding')
        self.assertEqual(manual_site.site_list[1].binding_model.width, 0.012)
        self.assertEqual(manual_site.site_list[2].binding_model.width, 0.01)
        self.assertEqual(manual_site.site_list[3].binding_model.width, 0.01)
        self.assertEqual(manual_site.site_list[4].binding_model, None)


