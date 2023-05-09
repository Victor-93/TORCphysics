import numpy as np


class Log:
    def __init__(self, size, frames, time, dt, structure, name, seed, site_list, initial_twist, initial_superhelical):
        self.metadata = []  # Will contain list of events
        self.size = size  # size of circuit
        self.frames = frames  # Total number of frames ran in simulation
        self.time = time  # Total simulation time
        self.dt = dt  # Timestep
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

    # TODO: Think about how to calculate the elongation rate: should be the time between (unbinding-binding)
    def site_overall_analysis(self):
        for i, site in enumerate(self.site_list):
            global_sum = False  # This variable is for enzymes that recognise bare DNA
            # And is used to update its global quantities.

            # This is for enzymes that bind bare DNA - Let's initialize it
            if 'DNA' in site.site_type and 'sum' == site.name:
                global_sum = True

            for event in self.metadata:
                if event.site.name != site.name and not global_sum:
                    continue
                if global_sum and event.site.site_type != site.site_type:
                    continue
                if event.event_type == 'binding_event':
                    self.total_binding_events[i] += 1
                if event.event_type == 'unbinding_event':
                    self.total_unbinding_events[i] += 1

        self.calculate_binding_rates()
        return

    def calculate_binding_rates(self):
        for i, site in enumerate(self.site_list):
            self.binding_rates[i] = self.total_binding_events[i]/self.time
        return

    # This one is the same as the binding rates.... Instead of unbinding should be elongation rates...
    # because then the unbinding rate would be the same as the binding rate
    def calculate_unbinding_rates(self):
        for i, site in enumerate(self.site_list):
            self.unbinding_rates[i] = self.total_unbinding_events[i]/self.time
        return

    def calculate_elongation_rates(self):
        pass

    # Writes information to a log file named self.name+.log
    def log_out(self):
        # TODO: Maybe you can clean a little bit your log (binding_rates), so you don't print EXT or extra stuff
        f = open(self.name + ".log", "w")
        f.write("TORCphysics log file \n")
        f.write("\n")
        f.write("General information: \n")
        f.write("---------------------------------------------------------------------------------------\n")
        f.write("Circuit name: " + self.name + "\n")
        f.write("Size: " + str(int(self.size)) + "\n")
        f.write("Structure: " + self.structure + "\n")
        f.write("Number of sites: " + str(int(self.n_sites)) + "\n")
        f.write("Number of frames: " + str(int(self.frames)) + "\n")
        f.write("Simulation time: " + str(self.time) + "\n")
        f.write("Simulation timestep: " + str(self.dt) + "\n")
        f.write("Seed: " + str(self.seed) + "\n")
        f.write(" \n")
        f.write("Overall simulation information \n")
        f.write("---------------------------------------------------------------------------------------\n")
        f.write("Initial twist: " + f'{self.initial_twist:.5f}\n')
        f.write("Final twist: " + f'{self.final_twist:.5f}' + "\n")
        f.write("Initial superhelical density: " + f'{self.initial_superhelical:.5f}' + "\n")
        f.write("Final superhelical density: " + f'{self.final_superhelical:.5f}' + "\n")
        f.write(" \n")
        f.write("Overall rates \n")
        f.write("---------------------------------------------------------------------------------------\n")
        # TODO: Better format in this one as well
        self.site_overall_analysis()
        c1 = 'site type'
        c2 = 'site name'
        c3 = '#binding events'
        c4 = '#unbinding events'
        c5 = 'binding rates'
        #c5 = 'Elongation rates'
        f.write(f'{c1:20} {c2:20} {c3:20} {c4:20} {c5:20} \n ')  # {c5:20} \n')
        for i, site in enumerate(self.site_list):
            f.write(f'{site.site_type:20} {site.name:20} {self.total_binding_events[i]:20} '
                    f'{self.total_unbinding_events[i]:20}  {self.binding_rates[i]:20}')
            # {self.elongation_rates[i]:20}')
            f.write(" \n")
        f.write(" \n")
        f.write("Events\n")
        f.write("---------------------------------------------------------------------------------------\n")
        # TODO: make a better format
        for event in self.metadata:
            # line = str(event.frame) + " " + str(event.time) + " " + event.event_type + " " + event.site.site_type +
            # \ " " + event.site.name + " " + event.enzyme.name + " " + str(event.twist) + " " + \ str(
            # event.superhelical) + " " + str(event.global_twist) + " " + str(event.global_superhelical) + "\n"
            # line
            # = 'frame ' + str(event.frame) + " time " + str(event.time) + " " + event.event_type + " at " + \
            # event.site.site_type + " " + event.site.name + " enzyme " + event.enzyme.name + " with twist " + \ str(
            # event.twist) + " superhelical " + str(event.superhelical) + " global twist " + \ str(
            # event.global_twist) + " global supercoiling " + str(event.global_superhelical) + "\n"
            line = 'frame ' + str(event.frame) + " time " + str(event.time) + " " + event.event_type + " at " + \
                   event.site.site_type + " " + event.site.name + " enzyme " + event.enzyme.name + " with twist " + \
                   f'{event.twist:.5f}' + " superhelical " + f'{event.superhelical:.5f}' + " global twist " + \
                   f'{event.global_twist:.5f}' + " global supercoiling " + f'{event.global_superhelical:.5f}' + "\n"
            f.write(line)
        f.close()


class Event:
    def __init__(self, time, frame, event_type, twist, superhelical, global_twist, global_superhelical,
                 site, enzyme, position):
        self.time = time
        self.frame = frame
        self.event_type = event_type
        self.twist = twist
        self.superhelical = superhelical
        self.global_twist = global_twist
        self.global_superhelical = global_superhelical
        self.site = site
        self.enzyme = enzyme
        self.position = position
