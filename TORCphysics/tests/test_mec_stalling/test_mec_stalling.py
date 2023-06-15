from unittest import TestCase
from TORCphysics import Circuit, add_enzyme, visualization


class TestMecstalling(TestCase):

    def test_torque_stall_Geng(self):
        circuit_filename = 'circuit.csv'
        sites_filename = 'sites.csv'
        enzymes_filename = 'enzymes.csv'
        environment_filename = 'environment.csv'
        output_prefix = 'output'
        frames = 500
        series = True
        continuation = False
        tm = 'stochastic'
        # tm = 'continuum'
        mm = 'torque_stall_Geng'
        dt = .5
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        my_site = [site for site in my_circuit.site_list if site.name == 'gene1'][0]
        my_env = [environment for environment in my_circuit.environmental_list if environment.name == 'RNAP'][0]

        RNAP1 = add_enzyme.AddEnzyme(20, e_type='RNAP', name='RNAP', site=my_site,
                                     position=my_site.start,
                                     size=my_env.size)
        RNAP2 = add_enzyme.AddEnzyme(200, e_type='RNAP', name='RNAP', site=my_site,
                                     position=my_site.start,
                                     size=my_env.size)
        custom_enzymes = [RNAP1, RNAP2]
        add_enzyme.run_simulation(my_circuit, custom_enzymes)

        visualization.create_animation_linear(my_circuit, my_circuit.sites_df,
                                              my_circuit.enzymes_df,
                                              my_circuit.name, '.gif', site_type='gene')
