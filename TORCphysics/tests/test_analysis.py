from unittest import TestCase
from TORCphysics import Circuit
from TORCphysics import analysis as an

class TestAnalysis(TestCase):

    def test_signal_by_type(self):
        circuit_filename = '../circuit.csv'
        sites_filename = '../sites.csv'
        enzymes_filename = '../enzymes.csv'
        environment_filename = '../environment.csv'
        output_prefix = 'output'
        frames = 1000
        series = True
        continuation = False
        tm = 'continuum'
        mm = 'uniform'
        dt = 1
        my_circuit = Circuit(circuit_filename, sites_filename, enzymes_filename, environment_filename,
                             output_prefix, frames, series, continuation, dt, tm, mm)
        my_circuit.run()
        my_signal = an.build_signal_by_type(my_circuit.sites_df, 'gene')

        self.assertGreater(my_circuit.get_num_enzymes(), 0, "Empty enzyme list")
        self.assertGreater(my_circuit.get_num_sites(), 0, "Empty enzyme list")
        self.assertGreater(my_circuit.get_num_environmentals(), 0, "Empty enzyme list")