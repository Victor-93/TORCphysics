from unittest import TestCase
from TORCphysics import Site, Enzyme
from TORCphysics import effect_model as em

class TestEffects(TestCase):

    def test_get_start_end_c(self):
        nbp = 5
        my_site = Site(s_type='test_site', name='u1', start=1, end=nbp, k_min=0, k_max=0, s_model_name=None, oparams=None)
        my_enzyme = Enzyme(e_type='test', name='YO', site=my_site, position=2, size=1, k_cat=0.0, k_off=0.0,
                           twist=0, superhelical=0)
        position_left, position_right = em.get_start_end_c(my_enzyme, my_enzyme, nbp)
        self.assertEqual(-2, position_left, "Incorrectly calculated the left position")
        self.assertEqual(7, position_right, "Incorrectly calculated the right position")

    def test_calculate_length(self):
        nbp = 2000
        my_site = Site(s_type='test_site', name='u1', start=1, end=nbp, k_min=0, k_max=0, s_model_name=None, oparams=None)
        enzyme1 = Enzyme(e_type='test', name='YO', site=my_site, position=200, size=50, k_cat=0.0, k_off=0.0,
                         twist=0, superhelical=0)
        enzyme2 = Enzyme(e_type='test', name='YO', site=my_site, position=1200, size=50, k_cat=0.0, k_off=0.0,
                         twist=0, superhelical=0)
        length = em.calculate_length(enzyme1, enzyme2)
        self.assertEqual(length, 950, "Incorrectly calculated the length")
