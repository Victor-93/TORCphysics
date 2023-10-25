from unittest import TestCase
from TORCphysics import models_workflow as mw
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as ubm
from TORCphysics import Site, Environment, Enzyme
import numpy as np


class TestWorkflows(TestCase):

    def test_doc(self):
        help(mw)

    def test_binding_workflow(self):
        # Let's create sites. We'll have 4 genes, two genes are type gene1 and the other two are
        # type gene2
        s0 = Site(site_type='EXT', name='EXT', start=1, end=5000, k_on=0.0)
        s1 = Site(site_type='gene1', name='gene1_1', start=100, end=500, k_on=0.0)
        s2 = Site(site_type='gene1', name='gene1_2', start=1100, end=1500, k_on=0.0)
        s3 = Site(site_type='gene2', name='gene2_1', start=2100, end=2500, k_on=0.0)
        s4 = Site(site_type='gene2', name='gene2_2', start=3100, end=3500, k_on=0.0)

        # And let's creat enzymes that act as boundaries
        extra_left = Enzyme(e_type='EXT', name='EXT_L', site=s0,
                            position=1, size=0, effective_size=0,
                            twist=0, superhelical=-0.06)
        extra_right = Enzyme(e_type='EXT', name='EXT_R', site=s0,
                             position=5000, size=0, effective_size=0,
                             twist=0, superhelical=-0.06)
        enzyme_list = [extra_left, extra_right]

        # And let's set our dt and rng
        dt = 1.0
        rng = np.random.default_rng(1993)  # my birth year

        # ---------------------------------------------------------------
        # Test 1. Environment with concentration = 0.0, returns empty list.
        # ---------------------------------------------------------------
        print('test1')
        s1.binding_model = bm.PoissonBinding(**{'k_on': 9999})
        s2.binding_model = bm.PoissonBinding(**{'k_on': 9999})
        s3.binding_model = bm.PoissonBinding(**{'k_on': 9999})
        s4.binding_model = bm.PoissonBinding(**{'k_on': 9999})
        site_list = [s0, s1, s2, s3, s4]

        # And let's create the environmentals. One of them recognise the first genes type 1 and the other type 2
        e1 = Environment(e_type='RNAP', name='test1', site_list=site_list, concentration=0.0,
                         size=60, effective_size=30, site_type='gene1')
        e2 = Environment(e_type='RNAP', name='test2', site_list=site_list, concentration=0.0,
                         size=60, effective_size=30, site_type='gene2')
        environment_list = [e1, e2]

        new_enzymes = mw.binding_workflow(enzyme_list=enzyme_list, environmental_list=environment_list, dt=dt, rng=rng)
        self.assertEqual(len(new_enzymes), 0)

        # ---------------------------------------------------------------
        # Test 2. Environment and site without binding model, returns empty list.
        # ---------------------------------------------------------------
        print('test2')
        s1.binding_model = None
        s2.binding_model = None
        s3.binding_model = None
        s4.binding_model = None
        # And let's create the environmentals. One of them recognise the first genes type 1 and the other type 2
        e1 = Environment(e_type='RNAP', name='test1', site_list=site_list, concentration=10.0,
                         size=60, effective_size=30, site_type='gene1')
        e2 = Environment(e_type='RNAP', name='test2', site_list=site_list, concentration=10.0,
                         size=60, effective_size=30, site_type='gene2')
        environment_list = [e1, e2]

        new_enzymes = mw.binding_workflow(enzyme_list=enzyme_list, environmental_list=environment_list, dt=dt, rng=rng)
        self.assertEqual(len(new_enzymes), 0)

        # ---------------------------------------------------------------
        # Test 3. Environment & sites but sites are global, returns empty list
        # ---------------------------------------------------------------
        print('test3')
        s1.binding_model = bm.PoissonBinding(**{'k_on': 1})
        s2.binding_model = bm.PoissonBinding(**{'k_on': 1})
        s3.binding_model = bm.PoissonBinding(**{'k_on': 1})
        s4.binding_model = bm.PoissonBinding(**{'k_on': 1})
        s1.global_site = True
        s2.global_site = True
        s3.global_site = True
        s4.global_site = True
        site_list = [s0, s1, s2, s3, s4]
        e1 = Environment(e_type='RNAP', name='test1', site_list=site_list, concentration=10.0,
                         size=60, effective_size=30, site_type='gene1')
        e2 = Environment(e_type='RNAP', name='test2', site_list=site_list, concentration=10.0,
                         size=60, effective_size=30, site_type='gene2')
        environment_list = [e1, e2]

        new_enzymes = mw.binding_workflow(enzyme_list=enzyme_list, environmental_list=environment_list, dt=dt, rng=rng)
        self.assertEqual(len(new_enzymes), 0)

        # ---------------------------------------------------------------
        # Test 4. Environment and sites with k_on=0, hence probability 0; returns empty list.
        # ---------------------------------------------------------------
        print('test4')
        s1.global_site = False
        s2.global_site = False
        s3.global_site = False
        s4.global_site = False
        s1.binding_model = bm.PoissonBinding(**{'k_on': 0})
        s2.binding_model = bm.PoissonBinding(**{'k_on': 0})
        s3.binding_model = bm.PoissonBinding(**{'k_on': 0})
        s4.binding_model = bm.PoissonBinding(**{'k_on': 0})
        site_list = [s0, s1, s2, s3, s4]
        e1 = Environment(e_type='RNAP', name='test1', site_list=site_list, concentration=10.0,
                         size=60, effective_size=30, site_type='gene1')
        e2 = Environment(e_type='RNAP', name='test2', site_list=site_list, concentration=10.0,
                         size=60, effective_size=30, site_type='gene2')
        environment_list = [e1, e2]
        new_enzymes = mw.binding_workflow(enzyme_list=enzyme_list, environmental_list=environment_list, dt=dt, rng=rng)
        self.assertEqual(len(new_enzymes), 0)

        # ---------------------------------------------------------------
        # Test 5. Environment and sites with k_on=1, hence probability 1; returns list with enzymes.
        # ---------------------------------------------------------------
        print('test5')
        s1.binding_model = bm.PoissonBinding(**{'k_on': 1})
        s2.binding_model = bm.PoissonBinding(**{'k_on': 1})
        s3.binding_model = bm.PoissonBinding(**{'k_on': 1})
        s4.binding_model = bm.PoissonBinding(**{'k_on': 1})
        site_list = [s0, s1, s2, s3, s4]
        e1 = Environment(e_type='RNAP', name='test1', site_list=site_list, concentration=10.0,
                         size=60, effective_size=30, site_type='gene1')
        e2 = Environment(e_type='RNAP', name='test2', site_list=site_list, concentration=10.0,
                         size=60, effective_size=30, site_type='gene2')
        environment_list = [e1, e2]
        new_enzymes = mw.binding_workflow(enzyme_list=enzyme_list, environmental_list=environment_list, dt=dt, rng=rng)
        self.assertGreater(len(new_enzymes), 0)  # With the rng, 2 enzymes bind. So don't change the seed!

        # ---------------------------------------------------------------
        # Test 6. Environment and sites with k_on=1, but sites occupied with enzymes, returns empty list.
        # ---------------------------------------------------------------
        print('test6')

        # And let's creat the enzymes that block the sites and add them.
        en1 = Enzyme(e_type='NAP', name='block1', site=s0, position=80, size=60, effective_size=60,
                     twist=0, superhelical=-0.06)
        en2 = Enzyme(e_type='NAP', name='block2', site=s0, position=1080, size=60, effective_size=60,
                     twist=0, superhelical=-0.06)
        en3 = Enzyme(e_type='NAP', name='block3', site=s0, position=2080, size=60, effective_size=60,
                     twist=0, superhelical=-0.06)
        en4 = Enzyme(e_type='NAP', name='block4', site=s0, position=3080, size=60, effective_size=60,
                     twist=0, superhelical=-0.06)
        enzyme_list = [extra_left, en1, en2, en3, en4, extra_right]

        new_enzymes = mw.binding_workflow(enzyme_list=enzyme_list, environmental_list=environment_list, dt=dt, rng=rng)
        self.assertEqual(len(new_enzymes), 0)

        # TODO: We need to test every part of the code.
        #  7.- Also check how enzymes are added.

    def test_effect_workflow(self):
        # TODO: Test cases:
        #  Test 1. Enzyme is EXT - returns empty list
        #  Test 2. Enzyme has no effect_model - returns empty list
        #  Test 3. Enzyme has effect_model - returns one Effect
        #  Test 4. List of environments with no effect_models - returns empty list
        #  Test 5. List of environments, but not continuum - returns empty list
        #  Test 6. List of continuum environments, but only EXT_ enzymes - returns empty list
        #  Test 7. List of M continuum environments, but with N enzymes - returns M*N Effects
        s=2

