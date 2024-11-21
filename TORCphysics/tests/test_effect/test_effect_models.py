from unittest import TestCase
from TORCphysics import Circuit
from TORCphysics import Site, Environment

circuit_circular_filename = '../test_circuit/circuit.csv'
circuit_linear_filename = '../test_circuit/circuit_linear.csv'
sites_filename = None#'sites.csv'  # Has one lacOP binding site
enzymes_filename = None#'../test_inputs/empty_enzymes.csv'
environment_filename = None#'environment.csv'  # Has specific binding enzymes and bare DNA binding. All with concentration=0
output_prefix = ''
frames = 2000
series = True
continuation = False
dt = 1
error = 0.0001

class TestRNAPStagesStall(TestCase):

    # So, we will define some sites (3) and run a simulation where RNAPs can bin but will not initiate, they may bind
    # unbind but they will not induce any supercoils.
    # Test: The tests consits in running simulations and checking the global superhelical density. It shouldn't change!

    def test_instant_twist_transfer(self):
        RNAP_environmental = Environment(e_type='RNAP', name='RNAP', site_list=[], concentration=1, size=30,
                                         effective_size=20, site_type='gene',
                                         effect_model_name='RNAPStagesStallv2',
                                         unbinding_model_name='RNAPStagesSimpleUnbindingv2')

        # These three genes will have elongation turned off
        site_gene1 = Site(site_type='gene', name='test_gene1', start=100, end=500, k_on=.1,binding_model_name='GaussianBinding')
        site_gene2 = Site(site_type='gene', name='test_gene2', start=600, end=800, k_on=0.1,binding_model_name='GaussianBinding')
        site_gene3 = Site(site_type='gene', name='test_gene3', start=1500, end=1000, k_on=0.1,binding_model_name='GaussianBinding')

        site_gene1.binding_model.k_ini = 0.0001 # Very small
        site_gene2.binding_model.k_ini = 0.0001 # Very small
        site_gene3.binding_model.k_ini = 0.0001 # Very small

        # And these two will be able to transcribe
        site_gene4 = Site(site_type='gene', name='test_gene4', start=4000, end=4200, k_on=0.1,binding_model_name='GaussianBinding')
        site_gene5 = Site(site_type='gene', name='test_gene5', start=5200, end=5000, k_on=0.1,binding_model_name='GaussianBinding')
        site_gene6 = Site(site_type='gene', name='test_gene5', start=5400, end=5600, k_on=0.1,binding_model_name='GaussianBinding')


        # Case 1: Circular and no elongation (but binding/unbinding enabled)
        # --------------------------------------------------------------------------------------------------------
        my_circuit = Circuit(circuit_circular_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt)

        # Let's append the sites
        site_list = [site_gene1, site_gene2, site_gene3]
        for site in site_list:
            my_circuit.site_list.append( site)
        RNAP_environmental.site_list = site_list
        my_circuit.environmental_list.append(RNAP_environmental)
        my_circuit.sort_lists()

        # Run simulation and get superhelical densities
        superhelical = my_circuit.run_return_global_supercoiling()

        # Test!
        diff = sum(superhelical)
        self.assertLessEqual(diff, error, "Superhelical changed")

        # Case 2: Linear and with three genes that can elongate!
        # TODO: For some reason, circular plasmid with elongating RNAPs fails! Are we doing something wrong here or
        #       is it a bug in the main code?
        # --------------------------------------------------------------------------------------------------------
        my_circuit = Circuit(circuit_linear_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt)

        # Let's append the sites
        site_list = [site_gene1, site_gene2, site_gene3, site_gene4, site_gene5, site_gene6]
        for site in site_list:
            my_circuit.site_list.append( site)
        RNAP_environmental.site_list = site_list
        my_circuit.environmental_list.append(RNAP_environmental)
        my_circuit.sort_lists()

        # Run simulation and get superhelical densities
        superhelical = my_circuit.run_return_global_supercoiling()

        # Test!
        diff = sum(superhelical)
        self.assertLessEqual(diff, error, "Superhelical changed")
