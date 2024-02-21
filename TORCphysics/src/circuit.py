import pandas as pd
import random
import numpy as np
import sys
from TORCphysics import (Site, SiteFactory, Enzyme, EnzymeFactory, Environment,
                         EnvironmentFactory, Event, Log, params, utils)
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm
from TORCphysics import models_workflow as mw


# TODO: Check which inputs are optional (this is not urgent right now).
class Circuit:

    def __init__(self, circuit_filename, sites_filename, enzymes_filename, environment_filename,
                 output_prefix, frames, series, continuation, dt, random_seed=random.randrange(sys.maxsize)):
        # I'll save the input filenames just in case
        self.circuit_filename = circuit_filename
        self.sites_filename = sites_filename
        self.enzymes_filename = enzymes_filename
        self.environment_filename = environment_filename
        self.output_prefix = output_prefix
        self.frames = frames
        self.frame = 0
        self.series = series
        self.continuation = continuation
        self.name = None
        self.structure = None
        self.circle = None
        self.size = 0
        self.twist = None
        self.superhelical = None
        self.sequence = None
        self.dt = dt  # time step
        self.read_csv()  # Here, it gets the name,structure, etc
        self.site_list = SiteFactory(filename=sites_filename).site_list
        self.enzyme_list = EnzymeFactory(filename=enzymes_filename, site_list=self.site_list).enzyme_list
        self.environmental_list = EnvironmentFactory(filename=environment_filename,
                                                     site_list=self.site_list).environment_list
        self.check_object_inputs()  # Checks that the site, environmental and enzyme lists are correct
        self.time = 0
        # create a time-based seed and save it, and initialize our random generator with this seed
        self.seed = random_seed  # random.randrange(sys.maxsize)
        self.rng = np.random.default_rng(self.seed)
        # This option indicates if you want to include the specific sites for non-specific DNA binding proteins
        self.write_nonspecific_sites = False  # TODO: add this option as input
        # TODO: Create list of warnings in the future, and remove duplicates to print at the end of simulation
        # self.warnings = []  # List with warnings

        # -----------------------------------------------
        # We add new DNA sites which is the ones that we will link topos binding
        # Note: In the future there might be cases in which new enzymes will be added to the environment, so maybe,
        # these new sites will need to be created
        # TODO: Maybe we can store it in define_bar_DNA_binding_sites, so this happens in general, for DNA binding
        #  enzymes.
        # Define bare DNA binding sites for bare DNA binding enzymes
        self.define_bare_DNA_binding_sites()

        # Sort list of enzymes and sites by position/start
        self.sort_lists()
        # Distribute twist/supercoiling
        self.add_fake_boundaries()
        self.sort_lists()

        self.update_global_twist()
        self.update_global_superhelical()

        # Let's initialize the log
        self.log = Log(self.size, self.frames, self.frames * self.dt, self.dt, self.structure,
                       self.name + '_' + output_prefix, self.seed,
                       self.site_list, self.twist, self.superhelical, self.write_nonspecific_sites)

        # Let's define the dictionaries that will become dataframes, in case the series option was selected
        self.enzymes_df = []
        self.enzymes_dict_list = []
        self.append_enzymes_to_dict()

        self.sites_df = []
        self.sites_dict_list = []
        self.sites_dict_list_aux = []  # this one is an auxiliary
        self.append_sites_to_dict_step1()
        self.append_sites_to_dict_step2([], [])

        self.environmental_df = []
        self.environmental_dict_list = []
        self.append_environmental_to_dict()

    # This reads the circuit csv and sorts out the twist and structure
    def read_csv(self):
        df = pd.read_csv(self.circuit_filename)
        sequence_file = df['sequence'][0]
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
        # And the sequence if given
        if (sequence_file == '' or sequence_file == 'None' or sequence_file == 'none' or sequence_file is None
                or sequence_file == np.nan):
            self.sequence = None
        else:
            self.sequence = utils.read_fasta(file_name=sequence_file)
            # self.read_fasta(sequence_file)
            self.size = len(self.sequence)

    # This one runs the simulation
    # TODO: Think a way of making your run() function run for additional number of frames. This will make the code more
    #  versatile and will allow you create experiments where you add stuff manually
    # TODO: It might be worth adding an action and time? Or not? So that maybe it could perform an action at a given
    #  time?
    def run(self):
        #  What I need to do for including more frames is modify the log as well, and all other places where
        #  self.frames is used...
        #  if n_frames is not None:
        #    frames_i = 1
        #    frames_f = n_frames + 1
        #  else:
        #    frames_i = 1
        #    frames_f = self.frames =
        #

        for frame in range(1, self.frames + 1):
            self.frame += 1
            self.time = frame * self.dt
            if self.series:
                self.append_sites_to_dict_step1()

            # BINDING
            # --------------------------------------------------------------
            # Apply binding model and get list of new enzymes
            new_enzyme_list = mw.binding_workflow(self.enzyme_list, self.environmental_list, self.dt, self.rng)
            # These new enzymes are lacking twist and superhelical, we need to fix them and actually add them
            # to the enzyme_list
            # But before, add the binding events to the log  (it's easier to do it first)
            #            self.add_binding_events_to_log(new_enzyme_list)
            self.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

            # EFFECT
            # --------------------------------------------------------------
            effects_list = mw.effect_workflow(self.enzyme_list, self.environmental_list, self.dt)
            self.apply_effects(effects_list)

            # UNBINDING
            # --------------------------------------------------------------
            drop_list_index, drop_list_enzyme = mw.unbinding_workflow(self.enzyme_list, self.dt, self.rng)
            self.drop_enzymes(drop_list_index)
            self.add_to_environment(drop_list_enzyme)
            #            self.add_unbinding_events_to_log(drop_list)

            # UPDATE GLOBALS
            # --------------------------------------------------------------
            self.update_global_twist()
            self.update_global_superhelical()

            # self.log.log_out()

            # Add to series df if the series option was selected (default=True)
            # --------------------------------------------------------------
            if self.series:
                self.append_enzymes_to_dict()
                self.append_sites_to_dict_step2(new_enzyme_list, drop_list_enzyme)
                self.append_environmental_to_dict()

        # Output the dataframes: (series)
        if self.series:
            self.enzymes_df = pd.DataFrame.from_dict(self.enzymes_dict_list)
            self.enzymes_df.to_csv(self.name + '_' + self.output_prefix + '_enzymes_df.csv', index=False, sep=',')
            self.sites_df = pd.DataFrame.from_dict(self.sites_dict_list)
            self.sites_df.to_csv(self.name + '_' + self.output_prefix + '_sites_df.csv', index=False, sep=',')
            self.environmental_df = pd.DataFrame.from_dict(self.environmental_dict_list)
            self.environmental_df.to_csv(self.name + '_' + self.output_prefix + '_environment_df.csv',
                                         index=False, sep=',')

        # Output the log of events
        self.log.final_twist = self.twist
        self.log.final_superhelical = self.superhelical
        self.log.log_out()

        # Output csvs
        self.enzyme_list_to_df().to_csv(self.name + '_' + self.output_prefix + '_enzymes' + '.csv',
                                        index=False, sep=',')
        self.site_list_to_df().to_csv(self.name + '_' + self.output_prefix + '_sites' + '.csv', index=False, sep=',')
        self.environmental_list_to_df().to_csv(self.name + '_' + self.output_prefix + '_environment' + '.csv',
                                               index=False, sep=',')

    # Sometimes, we might want to run a simulation but without producing any output file. Instead, we want
    # to run the simulation and return the actual dataframes so we can process them within the same script.
    # This function is particularly useful for that, as it returns the three types of dataframes: sites_df,
    # enzyme_df and environmental_df.
    def run_return_dfs(self):

        # run simulation
        for frame in range(1, self.frames + 1):
            self.frame += 1
            self.time = frame * self.dt
            if self.series:
                self.append_sites_to_dict_step1()

            # BINDING
            # --------------------------------------------------------------
            new_enzyme_list = mw.binding_workflow(self.enzyme_list, self.environmental_list, self.dt, self.rng)

            self.add_new_enzymes(new_enzyme_list)

            # EFFECT
            # --------------------------------------------------------------
            effects_list = mw.effect_workflow(self.enzyme_list, self.environmental_list, self.dt)
            self.apply_effects(effects_list)

            # UPDATE GLOBALS
            # --------------------------------------------------------------
            self.update_global_twist()
            self.update_global_superhelical()

            # UNBINDING
            # --------------------------------------------------------------
            drop_list_index, drop_list_enzyme = mw.unbinding_workflow(self.enzyme_list, self.dt, self.rng)
            self.drop_enzymes(drop_list_index)
            self.add_to_environment(drop_list_enzyme)

            # UPDATE GLOBALS
            # --------------------------------------------------------------
            self.update_global_twist()
            self.update_global_superhelical()

            # --------------------------------------------------------------
            if self.series:
                self.append_enzymes_to_dict()
                self.append_sites_to_dict_step2(new_enzyme_list, drop_list_enzyme)
                self.append_environmental_to_dict()

        # Output the dataframes: (series)
        if self.series:
            self.enzymes_df = pd.DataFrame.from_dict(self.enzymes_dict_list)
            self.sites_df = pd.DataFrame.from_dict(self.sites_dict_list)
            self.environmental_df = pd.DataFrame.from_dict(self.environmental_dict_list)

        return self.enzymes_df, self.sites_df, self.environmental_df

    # Sometimes we might be interested in the supercoiling global response, and not care about the specific interactions
    # This function performs a simulation where we do not save the bound/unbound enzymes, hence we do not produce
    # log, df or csv files. Only a numpy array -> my_supercoiling is returned.
    # This function is particularly useful when calibrating models with experiments, where it might be needed to run
    # hundreds of simulations and the calibration is performed according the global supercoiling responses to enzymes
    # such as topoisomerases.
    def run_return_global_supercoiling(self):
        my_supercoiling = np.zeros(self.frames + 1)
        my_supercoiling[0] = self.superhelical

        # run simulation
        for frame in range(1, self.frames + 1):
            # BINDING
            # --------------------------------------------------------------
            new_enzyme_list = mw.binding_workflow(self.enzyme_list, self.environmental_list, self.dt, self.rng)

            self.add_new_enzymes(new_enzyme_list)

            # EFFECT
            # --------------------------------------------------------------
            effects_list = mw.effect_workflow(self.enzyme_list, self.environmental_list, self.dt)
            self.apply_effects(effects_list)

            # UPDATE GLOBALS
            # --------------------------------------------------------------
            self.update_global_twist()
            self.update_global_superhelical()

            # UNBINDING
            # --------------------------------------------------------------
            drop_list_index, drop_list_enzyme = mw.unbinding_workflow(self.enzyme_list, self.dt, self.rng)
            self.drop_enzymes(drop_list_index)
            self.add_to_environment(drop_list_enzyme)

            # UPDATE GLOBALS
            # --------------------------------------------------------------
            self.update_global_twist()
            self.update_global_superhelical()
            my_supercoiling[frame] = self.superhelical
        return my_supercoiling

    # Returns list of enzymes in the form of dataframe. This function is with the intention of outputting the system
    def enzyme_list_to_df(self):
        enzyme_aux = []  # This will be a list of dicts
        for enzyme in self.enzyme_list:
            d = {'type': enzyme.enzyme_type, 'name': enzyme.name, 'site': enzyme.site.name, 'position': enzyme.position,
                 'size': enzyme.size, 'effective_size': enzyme.effective_size, 'twist': enzyme.twist,
                 'superhelical': enzyme.superhelical, 'effect_model': enzyme.effect_model_name,
                 'effect_oparams': enzyme.effect_oparams_file, 'unbinding_model': enzyme.unbinding_model_name,
                 'unbinding_oparams': enzyme.unbinding_oparams_file}
            enzyme_aux.append(d)
        my_df = pd.DataFrame.from_dict(enzyme_aux)
        return my_df

    # Returns environmental list in the form of dataframe. This function is with the intention of outputting the system
    def environmental_list_to_df(self):
        environmental_aux = []  # This will be a list of dicts
        for environmental in self.environmental_list:
            d = {'type': environmental.enzyme_type, 'name': environmental.name,
                 'site_type': environmental.site_type, 'concentration': environmental.concentration,
                 'size': environmental.size, 'effective_size': environmental.effective_size,
                 'binding_model': environmental.binding_model_name,
                 'binding_oparams': environmental.binding_oparams_file,
                 'effect_model': environmental.effect_model_name,
                 'effect_oparams': environmental.effect_oparams_file,
                 'unbinding_model': environmental.unbinding_model_name,
                 'unbinding_oparams': environmental.unbinding_oparams_file}
            environmental_aux.append(d)
        my_df = pd.DataFrame.from_dict(environmental_aux)
        return my_df

    # Returns list of sites in the form of dataframe. This function is with the intention of outputting the system
    def site_list_to_df(self):
        site_aux = []  # This will be a list of dicts
        for site in self.site_list:
            d = {'type': site.site_type, 'name': site.name, 'start': site.start, 'end': site.end, 'k_on': site.k_on,
                 'binding_model': site.binding_model_name, 'binding_oparams': site.binding_oparams_file}
            site_aux.append(d)
        my_df = pd.DataFrame.from_dict(site_aux)
        return my_df

    # Append new substances in the environment to the self.environmental_dict_list
    # These quantities are the frame, time, type of enzyme/substance, name and concentration
    def append_environmental_to_dict(self):
        for environmental in self.environmental_list:
            d = {'frame': self.frame, 'time': self.time, 'type': environmental.enzyme_type, 'name': environmental.name,
                 'concentration': environmental.concentration}
            self.environmental_dict_list.append(d)

    # Append new enzymes to the self.enzymes_dict_list
    # These quantities are at the end of the frame/time, where enzymes already bound/unbound and had an effect on the
    # circuit
    def append_enzymes_to_dict(self):
        for enzyme in self.enzyme_list:
            # if enzyme.enzyme_type == 'EXT':
            #    continue
            d = {'frame': self.frame, 'time': self.time, 'name': enzyme.name, 'site': enzyme.site.name,
                 'position': enzyme.position, 'twist': enzyme.twist, 'superhelical': enzyme.superhelical,
                 'global_twist': self.twist, 'global_superhelical': self.superhelical}
            self.enzymes_dict_list.append(d)

    # Append new useful data to the sites_dict_list.
    # This information corresponds to the events that happened during the current frame/time, where some enzymes
    # bound and unbound the DNA. The twist and superhelical density correspond to the ones before the binding happened,
    # so those two parameters do not correspond to the ones at the end of the time step
    # So this step1 function, it collects local twist and superhelical before binding. The step2 will count the
    # number of bound enzymes.
    def append_sites_to_dict_step1(self):
        self.sites_dict_list_aux.clear()  # Empty list
        d = {'frame': self.frame, 'time': self.time, 'type': 'circuit', 'name': self.name, 'twist': self.twist,
             'superhelical': self.superhelical, '#enzymes': 0,
             'binding': 0, 'unbinding': 0}
        self.sites_dict_list_aux.append(d)  # the first one is always the one corresponding to the circuit

        # Then collect local twist/supercoiling at each site before binding
        for site in self.site_list:
            if site.site_type == 'EXT':
                continue
            # skip non-specific binding proteins
            if not self.write_nonspecific_sites and site.name.isdigit() and 'DNA' in site.site_type:
                continue
            # This for global sites, e.g., for enzymes that bind bare DNA
            #  if 'DNA' in site.site_type and '_global' in site.name:
            if site.global_site:
                site_superhelical = self.superhelical
                site_twist = self.twist
            else:
                # TODO: Something is going on here with the EXT positions!
                enzyme_before = [enzyme for enzyme in self.enzyme_list if enzyme.position <= site.start]
                if len(enzyme_before) <= 0:
                    # TODO: Raises ValueError.
                    # TODO: The error that causes this, is that some enzymes go off the edges, find a way to prevent or
                    #  handle these situations
                    print('Some error in append_sites_dict_step1')
                    sys.exit()
                enzyme_before = [enzyme for enzyme in self.enzyme_list if enzyme.position <= site.start][-1]
                site_twist = enzyme_before.twist
                site_superhelical = enzyme_before.superhelical
            d = {'frame': self.frame, 'time': self.time, 'type': site.site_type, 'name': site.name, 'twist': site_twist,
                 'superhelical': site_superhelical, '#enzymes': 0, 'binding': 0, 'unbinding': 0}

            self.sites_dict_list_aux.append(d)

    # And the step2, where enzymes already bound/unbound. Here, it counts the number of enzymes that bound to each
    # site at the end of the frame. It also counts if during that frame, enzymes bound/unbound
    def append_sites_to_dict_step2(self, new_enzymes_list, drop_list_enzyme):
        # Let's modify the dictionary related to the whole circuit
        self.sites_dict_list_aux[0]['#enzymes'] = self.get_num_enzymes() - 2
        self.sites_dict_list_aux[0]['binding'] = len(new_enzymes_list)
        self.sites_dict_list_aux[0]['unbinding'] = len(drop_list_enzyme)
        self.sites_dict_list.append(self.sites_dict_list_aux[0])  # And add it to the big true list of dictionaries

        # Then collect local twist/supercoiling at each site before binding
        i = 0
        for site in self.site_list:
            if site.site_type == 'EXT':
                continue
            # skip non-specific binding proteins
            if not self.write_nonspecific_sites and site.name.isdigit() and 'DNA' in site.site_type:
                continue
            i = i + 1

            global_sum = False  # This variable is for enzymes that recognise bare DNA
            # And is used to update its global quantities.

            # This is for global sites, e.g., enzymes that bind bare DNA - Let's initialize it
            #  if 'DNA' in site.site_type and '_global' in site.name:
            if site.global_site:
                self.sites_dict_list_aux[i]['binding'] = 0
                self.sites_dict_list_aux[i]['unbinding'] = 0
                global_sum = True

            # Change binding to 1
            for new_enzyme in new_enzymes_list:
                if new_enzyme.site.name == site.name:
                    self.sites_dict_list_aux[i]['binding'] = 1
                # For globals in case of enzymes that bind bare DNA
                if global_sum and new_enzyme.site.site_type == site.site_type:
                    self.sites_dict_list_aux[i]['binding'] += 1

            # Change unbinding to 1
            for drop_enzyme in drop_list_enzyme:
                if drop_enzyme.site.name == site.name:
                    self.sites_dict_list_aux[i]['unbinding'] = 1
                # For globals in case of enzymes that bind bare DNA
                if global_sum and drop_enzyme.site.site_type == site.site_type:
                    self.sites_dict_list_aux[i]['unbinding'] += 1

            # This is mostly applied to genes, and checks how many enzymes are currently bound to that site
            self.sites_dict_list_aux[i]['#enzymes'] = \
                len([enzyme for enzyme in self.enzyme_list if enzyme.site.name == site.name])
            # And for the case of non-specific binding DNA proteins
            if global_sum:
                self.sites_dict_list_aux[i]['#enzymes'] = \
                    len([enzyme for enzyme in self.enzyme_list if enzyme.site.site_type == site.site_type])

            self.sites_dict_list.append(self.sites_dict_list_aux[i])  # And add it to the big true list of dictionaries

    # Calculates the global twist (just  sums the excess of twist)
    def update_global_twist(self):
        if self.circle:
            if self.get_num_enzymes() > 2:
                self.twist = sum(enzyme.twist for enzyme in self.enzyme_list) - self.enzyme_list[0].twist
            else:
                self.twist = sum(enzyme.twist for enzyme in self.enzyme_list)
        else:  # linear
            self.twist = sum(enzyme.twist for enzyme in self.enzyme_list)

    # And updates the global superhelical density
    # Important, assumes that global twist is updated
    def update_global_superhelical(self):
        if self.get_num_enzymes() > 2:
            self.superhelical = self.twist / (params.w0 * (self.size - sum(enzyme.size for enzyme in self.enzyme_list)))
        else:
            self.superhelical = self.twist / (params.w0 * self.size)

    def sort_lists(self):
        self.enzyme_list.sort(key=lambda x: x.position)
        self.site_list.sort(key=lambda x: x.start)

    def sort_site_list(self):
        self.site_list.sort(key=lambda x: x.start)

    def sort_enzyme_list(self):
        self.enzyme_list.sort(key=lambda x: x.position)

    #    def calculate_local_sites(self):
    #    #This function calculates the local sites with naked DNA.
    #    #If there are N enzymes, then there is N-1 local sites (or local domains).
    #        for i in range(self.get_num_enzymes()):
    #            self.site_list.append( Site())
    #        for enzyme in self.enzyme_list:
    #            print(0)

    def add_fake_boundaries(self):
        # I  need to add a fake site, so I can link the fake boundaries
        #    def __init__(self, site_type, name, start, end, k_on,
        #                 binding_model_name=None, binding_oparams_file=None,
        #                 binding_model=None, binding_model_oparams=None):
        # print(self.size)

        self.site_list.append(
            Site(site_type='EXT', name='EXT', start=1, end=float(self.size), k_on=0.0))

        # TODO: So the way you define continuations is with the fake boundaries? I should also include the local
        #  DNA sites
        if self.continuation:  # I need to fix this. It would be better if the output doesn't have EXT_L and EXT_R?
            a = 'EXT_L' in [x.name for x in self.enzyme_list]
            b = 'EXT_R' in [x.name for x in self.enzyme_list]
            if not (a and b):
                print('There is something wrong with the continuation file')
                print('Bye bye')
                sys.exit()
            else:
                print('Resuming simulation')

        else:  # If it is a new run
            if self.get_num_enzymes() > 0:  # This can only happen if there are objects bound to the DNA
                if self.circle:  # For circular DNA
                    position_left, position_right = utils.get_start_end_c(self.enzyme_list[0], self.enzyme_list[-1],
                                                                          self.size)

                else:  # For linear DNA
                    position_left = 0
                    position_right = self.size + 1

            else:  # If nothing is bound
                position_left = 0
                position_right = self.size + 1  # it is the same in this case for either linear or circular

            for enzyme in self.enzyme_list:  # Distribute supercoiling -
                enzyme.superhelical = self.superhelical

            # Let's treat the boundaries of our system as objects.
            # ----------------------------------------------------------------
            # Notice that we don't specify the end
            extra_left = Enzyme(e_type='EXT', name='EXT_L', site=self.site_match('EXT', 'EXT'),
                                position=float(position_left), size=0, effective_size=0,
                                twist=0, superhelical=self.superhelical)
            extra_right = Enzyme(e_type='EXT', name='EXT_R', site=self.site_match('EXT', 'EXT'),
                                 position=float(position_right), size=0, effective_size=0,
                                 twist=0, superhelical=self.superhelical)

            self.enzyme_list.append(extra_left)
            self.enzyme_list.append(extra_right)
            self.sort_lists()
            # And finally, update the twist
            self.update_twist()

            # WARNING!!!!
            # There could be a big mistake in case of linear structures that have a NAP in positions 1 or nbp

    # This functions updates twist in enzymes
    def update_twist(self):

        for i, enzyme in enumerate(self.enzyme_list[:-1]):
            enzyme.twist = utils.calculate_twist(enzyme, self.enzyme_list[i + 1])

    # Updates the supercoiling/superhelical in enzymes
    def update_supercoiling(self):
        for i, enzyme in enumerate(self.enzyme_list[:-1]):
            enzyme.superhelical = utils.calculate_supercoiling(enzyme, self.enzyme_list[i + 1])

    # Get number of enzymes
    def get_num_enzymes(self):
        return len(self.enzyme_list)

    # Gets number of environmentals
    def get_num_environmentals(self):
        return len(self.environmental_list)

    # Gets number of sites
    def get_num_sites(self):
        return len(self.site_list)

    # Matches labels with sites.
    def site_match(self, label_name, label_type):
        if label_name in [site.name for site in self.site_list] \
                and label_type in [site.site_type for site in self.site_list]:
            # TODO check if this works!
            for site in self.site_list:
                if site.name == label_name and site.site_type == label_type:
                    return site  # the first one?
        else:
            return None

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

    # Adds to the self.enzyme_list, the newly bound enzymes in new_enzyme_list
    # Also, creates the binding events and add them to the log. Notice that the twist and superhelical density are
    # the ones at the time of binding, before the effect and update
    # TODO: test this function!!!!
    def add_new_enzymes(self, new_enzyme_list):

        # Let's first sort the new list
        new_enzyme_list.sort(key=lambda x: x.position)

        #        print('before')
        #        print([enzyme.name for enzyme in self.enzyme_list])
        #        print([enzyme.twist for enzyme in self.enzyme_list])

        for new_enzyme in new_enzyme_list:

            # Get neighbour enzymes
            enzyme_before = [enzyme for enzyme in self.enzyme_list if enzyme.position <= new_enzyme.position][-1]
            enzyme_after = [enzyme for enzyme in self.enzyme_list if enzyme.position >= new_enzyme.position][0]

            # And quantities prior binding
            region_twist = enzyme_before.twist
            region_superhelical = enzyme_before.superhelical
            region_length = utils.calculate_length(enzyme_before, enzyme_after)

            # First, add the enzyme to the list and sort it
            self.enzyme_list.append(new_enzyme)
            self.sort_enzyme_list()

            # Before updating local parameters, create the new event and add it to log
            # --------------------------------------------------------------------------
            new_event = Event(self.time, self.frame, 'binding_event', enzyme_before.twist,
                              region_superhelical, self.twist, self.superhelical, new_enzyme.site, new_enzyme,
                              new_enzyme.position)
            # And add it to the log
            self.log.metadata.append(new_event)

            # Now we need to update the positions of the fake boundaries in circular DNA
            # --------------------------------------------------------------------------
            if self.circle:
                position_left, position_right = utils.get_start_end_c(self.enzyme_list[1], self.enzyme_list[-2],
                                                                      self.size)
                self.enzyme_list[0].position = position_left
                self.enzyme_list[-1].position = position_right

            # We are still missing the supercoiling density and the excess of twist...
            # We need to partition the twist, so it is conserved...
            # --------------------------------------------------------
            # These quantities are the sizes of the new local domains on the left and right of the new enzyme
            new_length_left = utils.calculate_length(enzyme_before, new_enzyme)
            new_length_right = utils.calculate_length(new_enzyme, enzyme_after)

            # now to calculate the new twists
            # NOTE that I don't partition using the supercoiling density because the region that is actually bound
            # is assumed to be relaxed by the enzyme. So the twist in the region increases because of the relaxed
            # bound region.
            # new_twist_left = region_twist * ((new_length_left + 0.5 * new_enzyme.size) / region_length)
            # new_twist_right = region_twist * ((new_length_right + 0.5 * new_enzyme.size) / region_length)
            new_twist_left = region_superhelical * region_length * new_length_left * params.w0 / (
                    new_length_left + new_length_right)
            new_twist_right = region_superhelical * region_length * new_length_right * params.w0 / (
                    new_length_left + new_length_right)
            # new_superhelical_left = new_twist_left / (params.w0*new_length_left)
            # new_superhelical_right = new_twist_right / (params.w0 * new_length_right)

            # update twists
            # ------------CIRCULAR DNA--------------------
            if self.circle:

                # There is no other domains besides the newly bound protein.
                if enzyme_before.name == 'EXT_L' and enzyme_after.name == 'EXT_R':
                    new_enzyme.twist = region_twist
                    # In this case, the twist of EXT_L and region_twist remain the same
                    # because it is a circular DNA with only one RNAP (no NAPs)

                # There is one EXT at the left
                elif enzyme_before.name == 'EXT_L' and enzyme_after.name != 'EXT_R':
                    # Check if this is how I can update a property in the enzymes - I think it does!
                    #                    enzyme_before.twist = new_twist_left
                    self.enzyme_list[0].twist = new_twist_left
                    self.enzyme_list[self.get_num_enzymes() - 2].twist = new_twist_left
                    # self.enzyme_list[-1].twist = new_twist_left
                    new_enzyme.twist = new_twist_right

                # There is one EXT at the right
                elif enzyme_before.name != 'EXT_L' and enzyme_before.name == 'EXT_R':
                    enzyme_before.twist = new_twist_left
                    self.enzyme_list[0] = new_twist_right
                    self.enzyme_list[-2] = new_twist_right
                    new_enzyme.twist = new_twist_right

                # In any other case where there's no neighbour boundaries
                else:
                    enzyme_before.twist = new_twist_left
                    new_enzyme.twist = new_twist_right

            # ------------LINEAR DNA--------------------
            else:
                enzyme_before.twist = new_twist_left
                new_enzyme.twist = new_twist_right

            # Now add the enzyme to the list, sort itself
            # self.enzyme_list.append(new_enzyme)
            # self.sort_enzyme_list()

            # And update supercoiling
            self.update_supercoiling()

    #        print('after')
    #        print([enzyme.name for enzyme in self.enzyme_list])
    #        print([enzyme.twist for enzyme in self.enzyme_list])

    # Drop enzymes specified in the drop_list. This list contains the indices in the self.enzyme_list that are unbinding
    # the DNA
    def drop_enzymes(self, drop_list):

        new_events = []  # List containing the new unbinding events

        for j, index in enumerate(drop_list):

            i = index - j  # This "i" is our true index. We subtract j because the indices (index) change by -1
            # everytime 1 enzyme is removed, hence we subtract -j

            # ------------CIRCULAR DNA--------------------
            if self.circle:  # As always, we need to be careful with the circular case

                # There is no other domains besides the newly bound protein. (O is the enzyme to be removed)
                # EXT_L........O..........EXT_R / The twist in O is passed to the left boundary (maybe a bit redundant?)
                if self.enzyme_list[i - 1].name == 'EXT_L' and self.enzyme_list[i + 1].name == 'EXT_R':
                    self.enzyme_list[0].twist = self.enzyme_list[i].twist  # Everything should have the same twist...?
                # There is one EXT at the left
                # EXT_L.............O........------------.....E......EXT_R / Twist in O has to be added to E,
                # and EXT_L becomes a mirrored version of E, so it has the same twist as E (which index is = N-2)
                elif self.enzyme_list[i - 1].name == 'EXT_L' and self.enzyme_list[i + 1].name != 'EXT_R':
                    #                    self.enzyme_list[self.get_num_enzymes() - 2].twist += self.enzyme_list[i].twist
                    #                    self.enzyme_list[0].twist = self.enzyme_list[self.get_num_enzymes() - 2].twist
                    self.enzyme_list[-2].twist += self.enzyme_list[i].twist
                    self.enzyme_list[0].twist = self.enzyme_list[-2].twist

                # ------.......E.......O.....---------- / Twist in O is added to E
                else:
                    self.enzyme_list[i - 1].twist += self.enzyme_list[i].twist
                    #                    self.enzyme_list[0].twist = self.enzyme_list[-2].twist
                    self.enzyme_list[0].twist = self.enzyme_list[i - 1].twist

            # ------------LINEAR DNA--------------------
            else:
                # ------.......E.......O.....---------- / Twist in O is added to E
                self.enzyme_list[i - 1].twist += self.enzyme_list[i].twist

            # Before removing the enzyme from the list, let's create the event and add it to the list of events
            # --------------------------------------------------------------------------
            new_event = Event(self.time, self.frame, 'unbinding_event', self.enzyme_list[i - 1].twist,
                              self.enzyme_list[i - 1].superhelical,
                              0, 0, self.enzyme_list[i].site, self.enzyme_list[i], self.enzyme_list[i].position)
            new_events.append(new_event)

            # Remove element of list
            # ------------------------------------------
            del self.enzyme_list[i]

        # Update fake boundaries positions if circular structure
        if self.enzyme_list[0].position > 0:
            # TODO: Check why EXT_L is changing its position
            this_error = self.enzyme_list[0]
            a = this_error
        if self.circle:
            if self.get_num_enzymes() > 2:
                self.enzyme_list[0].position, self.enzyme_list[-1].position = \
                    utils.get_start_end_c(self.enzyme_list[1], self.enzyme_list[-2], self.size)
            else:
                self.enzyme_list[0].position = 0
                self.enzyme_list[-1].position = self.size + 1

        # if self.enzyme_list[0].position > 0:
        #    print(0)
        self.sort_enzyme_list()
        self.update_supercoiling()

        # Now that the global supercoiling is updated, let's add the new unbinding events to the log
        for new_event in new_events:
            new_event.global_superhelical = self.superhelical
            new_event.global_twist = self.twist

            # And add it to the log
            self.log.metadata.append(new_event)

    # Apply effects in effects_list and realise output environment
    def apply_effects(self, effects_list):
        # And apply the effects for the specified enzymes in the effects_list
        for effect in effects_list:
            self.enzyme_list[effect.index].position += effect.position
            self.enzyme_list[effect.index].twist += effect.twist_right
            # In case we affect the boundary on the left - it affects the last (not fake) enzyme
            # because the fake boundaries don't exist and just reflect the first and last enzymes.
            if self.circle and effect.index == 1:
                self.enzyme_list[self.get_num_enzymes() - 2].twist += effect.twist_left
            if self.get_num_enzymes() > 2:  # In any other case just update the enzyme on the left
                self.enzyme_list[effect.index - 1].twist += effect.twist_left

        #            print(effect.index)
        #            print(0)
        # else:
        #    print('We have some issues in the effects, this should not be happening')
        #    sys.exit()

        # Now we need to update the positions of the fake boundaries in circular DNA
        # --------------------------------------------------------------------------
        if self.circle and self.get_num_enzymes() > 2:
            position_left, position_right = utils.get_start_end_c(self.enzyme_list[1],
                                                                  self.enzyme_list[self.get_num_enzymes() - 2],
                                                                  self.size)
            self.enzyme_list[0].position = position_left
            self.enzyme_list[-1].position = position_right
            self.enzyme_list[0].twist = self.enzyme_list[self.get_num_enzymes() - 2].twist

        # And update supercoiling - because twist was modified
        self.update_supercoiling()

    # Adds the output to the environment
    # TODO: there might be a better way to realise enzymes/substances to the environment.
    #  1.- Currently, only the concentration is summed. 2.- But will this still be the case if we add degradation?
    #  3.- And, will there be a more automatic way of defining these output to the environment?
    def add_to_environment(self, drop_list_enzymes):

        for enzyme in drop_list_enzymes:

            # TODO: Check this, it might notbe completly correct
            if enzyme.name == 'RNAP':
                size = abs(enzyme.site.start - enzyme.site.end + 1)
                output_environment = Environment(e_type='mRNA', name=enzyme.site.name, site_list=[], concentration=1,
                                                 size=size, effective_size=0, site_type=None)
            else:
                continue

            environment = [x for x in self.environmental_list if x.enzyme_type == output_environment.enzyme_type and
                           x.name == output_environment.name]
            if len(environment) > 0:
                environment[0].concentration += 1  # Maybe is not concentration in this case...
            else:
                self.environmental_list.append(output_environment)

    # This function define topoisomerase binding sites when using the stochastic binding model.
    # The way it works is that, it goes through the topoisomerases in the environment, then checks the empty space
    # between the enzymes O_____________________________O, then divides these empty spaces into binding sites in which
    # the topoisomerases would fit...
    def define_topoisomerase_binding_sites(self):
        topo_list = [environment for environment in self.environmental_list
                     if environment.enzyme_type == 'topo' or environment.enzyme_type == 'topoisomerase']
        for topo in topo_list:
            s = 0
            for i, enzyme in enumerate(self.enzyme_list[:-1]):
                next_enzyme = self.enzyme_list[i + 1]
                length = utils.calculate_length(enzyme, next_enzyme)
                n_sites = int(length / topo.size)
                for n in range(n_sites):  # The 1+n is to leave some space 1 empty space between enzymes
                    start = enzyme.position + enzyme.size + topo.size * n + (1 + n)
                    end = 1 + enzyme.position + enzyme.size + topo.size * (1 + n)
                    if end > next_enzyme.position:  # Little break to avoid the overlapping of enzymes
                        continue
                    topo_site = Site(s_type='DNA_' + topo.name, name=str(s), start=start, end=end, k_min=0, k_max=0,
                                     s_model_name='stochastic_' + topo.name, oparams=None)
                    self.site_list.append(topo_site)

                    s = s + 1
                    topo.site_list.append(topo_site)
        self.sort_site_list()
        return

    # TODO: Comment, test and sort.
    # This function defines the binding sites of enzymes that recognize bare DNA, that means just DNA.
    # It partitions the DNA in N binding sites of size enzyme.size
    def define_bare_DNA_binding_sites(self):

        # Let's create the global sites.
        # -----------------------------------------------

        # These global sites count how many times enzymes bound to the DNA molecule in general, but these sites
        # are not actually bound for any type of environmental, so they don't have binding/effect/unbinding models.
        for environmental in self.environmental_list:
            if 'DNA' in environmental.site_type:
                t_site = Site(site_type='DNA_' + environmental.name, name='DNA_' + environmental.name,
                              start=1, end=float(self.size), k_on=0.0, global_site=True)
                self.site_list.append(t_site)

        # Let's create the local sites.
        # -----------------------------------------------
        for environmental in self.environmental_list:
            if environmental.binding_model is None:  # There is no point in defining local sites if the environmentals
                # don't have a binding model
                continue
            if 'DNA' in environmental.site_type:
                n_sites = int(self.size / environmental.size)
                s = 0
                for n in range(n_sites):
                    start = 1 + environmental.size * n
                    end = environmental.size * (1 + n)
                    if end > self.size:  # Little break to avoid making it bigger than the actual plasmid
                        continue
                    local_site = Site(site_type='DNA_' + environmental.name, name=str(s), start=start, end=end,
                                      k_on=environmental.binding_model.k_on, binding_model=environmental.binding_model)
                    self.site_list.append(local_site)
                    environmental.site_list.append(local_site)

                    s = s + 1

                # The next line makes the environmental recognize the specific binding site
                environmental.site_type = 'DNA_' + environmental.name

        self.sort_site_list()
        return

    def check_object_inputs(self):

        """
         Checks that the input objects (sites, environments, enzymes) have the correct parameters or if they were not
         provided, it assigns them to the default models if necessary
         """

        # TODO: The sites, environmentals and enzymes check their inputs on their own. Here we should focus on checking
        #  the inputs of the circle class
        # TODO: Also check that size of bound enzymes, or sites are within the size of the circuit.
        # Sites - The Site class is capable of checking this on its own.
        # ================================
        # for site in self.site_list:
        #    if site.site_model == 'maxmin' and site.oparams is None:
        #        site.site_model = 'sam'
        #        print('For site ', site.name, 'maxmin model selected but no parameters provided. '
        #                                      'Switching to sam model')

        # Environment
        # ================================
        # for environment in self.environmental_list:

        # For the binding model
        # -------------------------------------------------------------------
        #   if environment.binding_oparams is None:#

        # In case of topoisomerase enzyme and stochastic model
        #       if self.topoisomerase_model == 'stochastic' and 'topo' in environment.enzyme_type:
        #           environment.binding_model = 'recognition'  # The default recognition model in case of stochastic
        #           if 'topoI' in environment.name:
        #               environment.binding_oparams = {'width': params.topo_b_w, 'threshold': params.topo_b_t}
        #           elif 'gyrase' in environment.name:
        #               environment.binding_oparams = {'width': params.gyra_b_w, 'threshold': params.gyra_b_t}

        # In case it is a continuum model, then there's no binding
        #  if self.topoisomerase_model == 'continuum' and 'topo' in environment.enzyme_type:
        #      environment.binding_model = None
        #      environment.binding_oparams = None

        # For the effect model
        # -------------------------------------------------------------------
        #   if environment.effect_oparams is None:
        # In case of topo
        #       if 'topo' in environment.enzyme_type:
        #           if 'topoI' in environment.name:
        #               environment.effect_oparams = {'k_cat': params.topo_k, 'width': params.topo_e_w,
        #                                             'threshold': params.topo_e_t}
        #           elif 'gyrase' in environment.name:
        #               environment.effect_oparams = {'k_cat': params.gyra_k, 'width': params.gyra_e_w,
        #                                             'threshold': params.gyra_e_t}
        # In case of RNAP
        #       elif 'RNAP' in environment.enzyme_type:
        #           environment.effect_oparams = {'velocity': params.v0, 'gamma': params.gamma}

        # Enzyme
        # ================================
        # for enzyme in self.enzyme_list:
        #    if enzyme.effect_oparams is None:
        # In case of topo
        #        if 'topo' in environment.enzyme_type:
        #            if self.topoisomerase_model == 'continuum':
        #                enzyme.effect_model = 'continuum'
        #            if 'topoI' in environment.name:
        #                environment.effect_oparams = {'k_cat': params.topo_k, 'width': params.topo_e_w,
        #                                              'threshold': params.topo_e_t}
        #            elif 'gyrase' in environment.name:
        #                environment.effect_oparams = {'k_cat': params.gyra_k, 'width': params.gyra_e_w,
        #                                              'threshold': params.gyra_e_t}
        # In case of RNAP
        #        elif 'RNAP' in environment.enzyme_type:
        #            environment.effect_oparams = {'velocity': params.v0, 'gamma': params.gamma}
