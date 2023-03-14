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
        frames = 200
        series = True
        continuation = False
        dt = 1
        tm = 'continuum'
        mm = 'uniform'
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        for site in my_circuit.site_list:  # I'll increase the rates
            # site.k_min = site.k_min * 100
            site.k_min = site.k_min * 10
        # Let's make the rates
        my_circuit.run()
        print(0)
