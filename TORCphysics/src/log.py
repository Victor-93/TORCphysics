import numpy as np


class Log:
    def __init__(self, size, frames, time, structure, name, seed, site_list, initial_twist, initial_superhelical):
        self.metadata = []  # Will contain list of events
        self.size = size  # size of circuit
        self.frames = frames  # Total number of frames ran in simulation
        self.time = time  # Total simulation time
        self.structure = structure  # Structure - linear or circular
        self.name = name  # name of circuit/experiment
        self.seed = seed
        self.initial_superhelical = initial_superhelical  # Global superhelical density
        self.initial_twist = initial_twist  # Global excess of twist
        self.site_list = site_list
        self.n_sites = len(site_list)  # number of sites

        # We still don't know these two last quantities
        self.final_superhelical = 0  # Global superhelical density
        self.final_twist = 0  # Global excess of twist

        # Useful information per site - Calculated with the metadata
        self.total_unbinding_events = np.zeros(self.n_sites)
        self.total_binding_events = np.zeros(self.n_sites)
        self.binding_rates = np.zeros(self.n_sites)
        self.unbinding_rates = np.zeros(self.n_sites)
        self.elongation_rates = np.zeros(self.n_sites)

        # Useful information per enzymes
        # total number of enzymes bound by type?

    def calculate_unbinding_events(self):
        pass

    def calculate_binding_events(self):
        pass

    def calculate_binding_rates(self):
        pass

    def calculate_unbinding_rates(self):
        pass

    def calculate_elongation_rates(self):
        pass

    # Writes information to a log file named self.name+.log
    def log_out(self):
        # TODO: Maybe you can clean a little bit your log (binding_rates), so you don't print EXT or extra stuff
        f = open(self.name+".log", "w")
        f.write("TORCphysics log file \n")
        f.write("\n")
        f.write("General information: \n")
        f.write("---------------------------------------------------------------------------------------\n")
        f.write("Circuit name: " + self.name +"\n")
        f.write("Size: \n")
        f.write("Structure: \n")
        f.write("Number of sites: \n")
        f.write("Number of frames: \n")
        f.write("Simulation time: \n")
        f.write("Seed:\n")
        f.write(" \n")
        f.write("Events\n")
        f.write("---------------------------------------------------------------------------------------\n")
#        '{:-9} YES votes  {:2.2%}'.format(yes_votes, percentage)

        for event in self.metadata:
#            line = str(event.frame) + " " + str(event.time) + " " + event.event_type + " " + event.site.site_type + \
 #                  " " + event.site.name + " " + event.enzyme.name + " " + str(event.twist) + " " + \
 #                  str(event.superhelical) + " " + str(event.global_twist) + " " + str(event.global_superhelical) + "\n"
            line = 'frame ' + str(event.frame) + " time " + str(event.time) + " " + event.event_type + " at " + \
                   event.site.site_type + " " + event.site.name + " enzyme " + event.enzyme.name + " with twist " + \
                   str(event.twist) + " superhelical " + str(event.superhelical) + " global twist " + \
                   str(event.global_twist) + " global supercoiling " + str(event.global_superhelical) + "\n"
            f.write(line)
        f.close()
