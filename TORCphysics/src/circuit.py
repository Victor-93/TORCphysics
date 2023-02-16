import pandas as pd
import random
import sys
from TORCphysics import SiteFactory, EnzymeFactory, EnvironmentFactory


class Circuit:

    def __init__(self, circuit_filename, sites_filename, enzymes_filename, environment_filename,
                 output_prefix, frames, series, continuation):
        # I'll save the input filenames just in case
        self.circuit_filename = circuit_filename
        self.sites_filename = sites_filename
        self.enzymes_filename = enzymes_filename
        self.environment_filename = environment_filename
        self.output_prefix = output_prefix
        self.frames = frames
        self.series = series
        self.continuation = continuation
        self.name = None
        self.structure = None
        self.circle = None
        self.size = None
        self.twist = None
        self.superhelical = None
        self.read_csv()  # Here, it gets the name,structure, etc
        self.sites = SiteFactory(sites_filename).site_list
        self.enzymes = EnzymeFactory(enzymes_filename, self.sites).enzyme_list
        self.environmentals = EnvironmentFactory(environment_filename, self.sites).environment_list
        self.time = 0
        # create a time-based seed and save it
        self.seed = random.randrange(sys.maxsize)

    # This reads the circuit csv and sorts out the twist and structure
    def read_csv(self):
        df = pd.read_csv(self.circuit_filename)
        self.name = df['name'][0]
        self.structure = df['structure'][0]
        if self.structure == 'circle' or self.structure == 'circular' or self.structure == 'close':
            self.circle = True
        elif self.structure == 'linear' or self.structure == 'lineal' or self.structure == 'open':
            self.circle = False
        else:
            self.circle = False
        self.size = df['size'][0]
        self.twist = df['twist'][0]
        self.superhelical = df['superhelical'][0]
        # TODO: Update twist - even if it's provided, it has to match the supercoiling density
        # (Maybe it shouldn't be provided)

    # Get number of enzymes
    def get_num_enzymes(self):
        return len(self.enzymes)

    # Gets number of environmentals
    def get_num_environmentals(self):
        return len(self.environmentals)

    # Gets number of sites
    def get_num_sites(self):
        return len(self.sites)

    # Prints general information
    def print_general_information(self):
        print("Running simulation")
        if self.circle:
            print("Circular structure")
        else:
            print("Linear structure")
        print("Running {0} frames on system composed of {1} bp".format(self.frames, self.size))
        print("Initial supercoiling density: {0}".format(self.superhelical))
        print("Initial twist: {0}".format(self.twist))
        print("Number of sites: {0}".format(self.get_num_sites()))
        print("Initial number of bound enzymes: {0}".format(self.get_num_enzymes()))
        print("Number of environmentals: {0}".format(self.get_num_environmentals()))
        print("Random seed: {0}".format(self.seed))

