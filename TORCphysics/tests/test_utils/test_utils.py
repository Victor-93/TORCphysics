from unittest import TestCase
from TORCphysics import utils, Enzyme, Site, Environment
import numpy as np


# TODO: Test remaining functions of utils.
class TestUtilsFunctions(TestCase):

    def test_instant_twist_transfer(self):


        # Here, we have a system in which an enzyme E is sorrouned by left L and right barrriers as: L__E____R
        # Enzyme E does not form a barrier and instantly lets the twist pass through.

        # Let's define objects we need
        # ------------------------------------------------
        # This a test site to link enzymes
        s0 = Site(site_type='EXT', name='EXT', start=1, end=5000, k_on=0.0)

        E_right = Enzyme(e_type='RNAP', name='RIGHT', site=s0, size=100, effective_size=50, position=2000,
                         twist=0.0, superhelical=0.0)

        # (Null). Both sides have the same superhelical density initially, so the change in twist is 0 on
        #         both sides
        # Total twist = 0 here
        E_left = Enzyme(e_type='RNAP', name='RIGHT', site=s0, size=100, effective_size=50, position=1000,
                         twist=0.0, superhelical=0.0)
        enzyme = Enzyme(e_type='RNAP', name='test1', site=s0, size=100, effective_size=50, position=1500,
                         twist=0.0, superhelical=0.0)
        enzyme_list = [E_left, enzyme, E_right]
        dtwist_left, dtwist_right = utils.instant_twist_transfer(enzyme,enzyme_list)
        self.assertEqual(dtwist_left, 0.0, 'Error in instant twist transfer')
        self.assertEqual(dtwist_right, 0.0, 'Error in instant twist transfer')

        # Same superhelical but not 0!
        E_left = Enzyme(e_type='RNAP', name='RIGHT', site=s0, size=100, effective_size=50, position=1000,
                         twist=-100.0, superhelical=0.0)
        enzyme = Enzyme(e_type='RNAP', name='test1', site=s0, size=100, effective_size=50, position=1500,
                         twist=-100.0, superhelical=0.0)
        enzyme_list = [E_left, enzyme, E_right]
        dtwist_left, dtwist_right = utils.instant_twist_transfer(enzyme,enzyme_list)
        self.assertEqual(dtwist_left, 0.0, 'Error in instant twist transfer')
        self.assertEqual(dtwist_right, 0.0, 'Error in instant twist transfer')

        # Negative superhelical density
        E_left = Enzyme(e_type='RNAP', name='RIGHT', site=s0, size=100, effective_size=50, position=1000,
                         twist=-100.0, superhelical=0.0)
        enzyme = Enzyme(e_type='RNAP', name='test1', site=s0, size=100, effective_size=50, position=1500,
                         twist=-1500.0, superhelical=0.0)
        enzyme_list = [E_left, enzyme, E_right]
        total_twist = E_left.twist + enzyme.twist
        dtwist_left, dtwist_right = utils.instant_twist_transfer(enzyme,enzyme_list)
        diff = abs(dtwist_left + E_left.twist + dtwist_right + enzyme.twist - total_twist)
        self.assertEqual(diff, 0.0, 'Error in instant twist transfer')

        # Negative and positive
        E_left = Enzyme(e_type='RNAP', name='RIGHT', site=s0, size=100, effective_size=50, position=1000,
                         twist=-200.0, superhelical=0.0)
        enzyme = Enzyme(e_type='RNAP', name='test1', site=s0, size=100, effective_size=50, position=1500,
                         twist=1000.0, superhelical=0.0)
        enzyme_list = [E_left, enzyme, E_right]
        total_twist = E_left.twist + enzyme.twist
        dtwist_left, dtwist_right = utils.instant_twist_transfer(enzyme,enzyme_list)
        diff = abs(dtwist_left + E_left.twist + dtwist_right + enzyme.twist - total_twist)
        self.assertEqual(diff, 0.0, 'Error in instant twist transfer')

        # Positive superhelical density
        E_left = Enzyme(e_type='RNAP', name='RIGHT', site=s0, size=100, effective_size=50, position=1000,
                         twist=2000.0, superhelical=0.0)
        enzyme = Enzyme(e_type='RNAP', name='test1', site=s0, size=100, effective_size=50, position=1500,
                         twist=1000.0, superhelical=0.0)
        enzyme_list = [E_left, enzyme, E_right]
        total_twist = E_left.twist + enzyme.twist
        dtwist_left, dtwist_right = utils.instant_twist_transfer(enzyme,enzyme_list)
        diff = abs(dtwist_left + E_left.twist + dtwist_right + enzyme.twist - total_twist)
        self.assertEqual(diff, 0.0, 'Error in instant twist transfer')

class TestSiteAvailability(TestCase):

    def test_site_availability(self):
        # If site is gene, RNAPs bind behind the site start:
        #  __RNAP_START____END, for direction = 1
        #  __END_______START_RNAP____, for direction = -1

        # If site is not a gene, hence, does not have direction, Enzymes (like NAPs) bind on the start site.
        #  _____START-NAP____END____ or ____END____START-NAP____ - These sites do not need an END anyway

        #  M = Molecule
        #  Test 1. direction = 1 _____M____START_____M_____END : Enzymes do not block the site
        #  Test 2. direction = 1 _____M____M-START-M____M_____END : One Enzyme blocks the site
        #  Test 3. direction = -1 _____M____END____M_____START___ : Enzymes do not block the site
        #  Test 4. direction = -1 _____M____END____M_____M-START-M___ : One Enzyme blocks the site
        #  Test 5. direction = 0 _____M____START_____M_____END : Enzymes do not block the site
        #  Test 6. direction = 0 _____M____M-START-M____M_____END : One Enzyme blocks the site

        # Let's create boundaries
        s0 = Site(site_type='EXT', name='EXT', start=1, end=5000, k_on=0.0)

        # Let's creat enzymes that act as boundaries
        extra_left = Enzyme(e_type='EXT', name='EXT_L', site=s0,
                            position=1, size=0, effective_size=0,
                            twist=0, superhelical=-0.06)
        extra_right = Enzyme(e_type='EXT', name='EXT_R', site=s0,
                             position=5000, size=0, effective_size=0,
                             twist=0, superhelical=-0.06)

        my_site = Site(site_type='gene1', name='test', start=100, end=500, k_on=0.0)
        # e_bind = Enzyme(e_type='RNAP', name='RNAP', site=my_site, position=80, size=60, effective_size=30,
        #                twist=0, superhelical=-0.06)

        e_bind = Environment(e_type='RNAP', name='test1', site_list=[my_site], concentration=10.0,
                             size=60, effective_size=30, site_type='gene1')

        # Test 1. direction = 1 _____M____START_____M_____END : Enzymes do not block the site
        enzyme_list = [extra_left, extra_right]
        available = utils.check_site_availability(site=my_site, enzyme_list=enzyme_list, environmental=e_bind)
        self.assertTrue(available)

        # Test 2. direction = 1 _____M____M-START-M____M_____END : One Enzyme blocks the site
        e_block = Enzyme(e_type='NAP', name='block1', site=my_site, position=80, size=60, effective_size=30,
                         twist=0, superhelical=-0.06)
        enzyme_list = [extra_left, e_block, extra_right]
        available = utils.check_site_availability(site=my_site, enzyme_list=enzyme_list, environmental=e_bind)
        self.assertFalse(available)

        #  Test 3. direction = -1 _____M____END____M_____START___ : Enzymes do not block the site
        my_site = Site(site_type='gene1', name='test', start=500, end=100, k_on=0.0)
        enzyme_list = [extra_left, extra_right]
        available = utils.check_site_availability(site=my_site, enzyme_list=enzyme_list, environmental=e_bind)
        self.assertTrue(available)

        #  Test 4. direction = -1 _____M____END____M_____M-START-M___ : One Enzyme blocks the site
        e_block = Enzyme(e_type='NAP', name='block1', site=my_site, position=505, size=60, effective_size=30,
                         twist=0, superhelical=-0.06)
        enzyme_list = [extra_left, e_block, extra_right]
        available = utils.check_site_availability(site=my_site, enzyme_list=enzyme_list, environmental=e_bind)
        self.assertFalse(available)

        #  Test 5. direction = 0 _____M____START_____M_____END : Enzymes do not block the site
        my_site = Site(site_type='IHF', name='test', start=500, end=0, k_on=0.0)
        e_bind = Environment(e_type='NAP', name='test1', site_list=[my_site], concentration=10.0,
                             size=100, effective_size=50, site_type='gene1')
        enzyme_list = [extra_left, extra_right]
        available = utils.check_site_availability(site=my_site, enzyme_list=enzyme_list, environmental=e_bind)
        self.assertTrue(available)

        #  Test 6. direction = 0 _____M____M-START-M____M_____END : One Enzyme blocks the site
        e_block = Enzyme(e_type='NAP', name='block1', site=my_site, position=470, size=60, effective_size=30,
                         twist=0, superhelical=-0.06)
        enzyme_list = [extra_left, e_block, extra_right]
        available = utils.check_site_availability(site=my_site, enzyme_list=enzyme_list, environmental=e_bind)
        self.assertFalse(available)
