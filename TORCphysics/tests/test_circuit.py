from unittest import TestCase
from TORCphysics import Circuit


class TestCircuit(TestCase):

    # For this test, Circuit should be able to process these non-empty input files (which hopefully are correct).
    # Then, it performs three tests: it counts the number of enzymes, sites and environmentals.
    def test_Circuit_init(self):
        circuit_filename = '../circuit.csv'
        sites_filename = '../sites.csv'
        enzymes_filename = '../enzymes.csv'
        environment_filename = '../environment.csv'
        output_prefix = 'output'
        frames = 5
        series = True
        continuation = False
        tm = 'continuum'
        mm = 'uniform'
        dt = 1
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        self.assertGreater(my_circuit.get_num_enzymes(), 0, "Empty enzyme list")
        self.assertGreater(my_circuit.get_num_sites(), 0, "Empty enzyme list")
        self.assertGreater(my_circuit.get_num_environmentals(), 0, "Empty enzyme list")

    # TODO: This can be several tests in one? 1.- Prove that statistical binding works, 2.- That the effects take place
    #  3.- And that the proteins unbind at the end?
    def test_run(self):
        circuit_filename = '../circuit.csv'
        sites_filename = '../sites.csv'
        enzymes_filename = '../enzymes.csv'
        environment_filename = '../environment.csv'
        output_prefix = 'output'
        frames = 1500
        series = True
        continuation = False
        dt = 0.5
        tm = 'continuum'
        mm = 'uniform'
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
#        for site in my_circuit.site_list:  # I'll increase the rates
#             site.k_min = site.k_min * 100
            #site.k_min = site.k_min * 5
        # Let's make the rates
        my_circuit.run()
        print(0)

    def test_run2(self):
        circuit_filename = '../circuit.csv'
        sites_filename = '../sites.csv'
        enzymes_filename = '../enzymes.csv'
        environment_filename = '../environment.csv'
        output_prefix = 'output'
        frames = 1500
        series = True
        continuation = False
        dt = 0.5
        tm = 'continuum'
        mm = 'uniform'
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        my_circuit.environmental_list[0].concentration = 0
        my_circuit.environmental_list[1].concentration = 0
        print(my_circuit.environmental_list)
        for site in my_circuit.site_list:  # I'll increase the rates
            # site.k_min = site.k_min * 100
            site.k_min = site.k_min * 0
        print(my_circuit.site_list[2].name)
        my_circuit.site_list[2].k_min = 0.01
        # Let's make the rates
        my_circuit.run()
        print(0)

    # TODO: test if enzymes bind correctly just before the start sites
    # TODO: test where one molecule binds, and the global superhelical must remain the same. create the system with 1
    #  gene, where only one RNAP binds, and unbinds. the global should remain the same.
    #  Maybe create a custom case using the functions in the circuit.

#    def test_add_1_enzyme(self):
