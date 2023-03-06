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
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, 1)
        self.assertGreater(my_circuit.get_num_enzymes(), 0, "Empty enzyme list")
        self.assertGreater(my_circuit.get_num_sites(), 0, "Empty enzyme list")
        self.assertGreater(my_circuit.get_num_environmentals(), 0, "Empty enzyme list")
