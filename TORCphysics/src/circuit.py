import pandas as pd
import random
import numpy as np
import sys
from TORCphysics import Site, SiteFactory, Enzyme, EnzymeFactory, Environment, EnvironmentFactory, Event, Log, params
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm


# TODO: Check how you determine enzymes bind, do they bind where with the left end on the site.start?
class Circuit:

    def __init__(self, circuit_filename, sites_filename, enzymes_filename, environment_filename,
                 output_prefix, frames, series, continuation, dt, topoisomerase_model, mechanical_model):
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
        self.dt = dt  # time step
        self.topoisomerase_model = topoisomerase_model
        self.mechanical_model = mechanical_model
        self.read_csv()  # Here, it gets the name,structure, etc
        self.site_list = SiteFactory(sites_filename).site_list
        self.enzyme_list = EnzymeFactory(enzymes_filename, self.site_list).enzyme_list
        self.environmental_list = EnvironmentFactory(environment_filename, self.site_list).environment_list
        self.check_object_inputs()  # Checks that the site, environmental and enzyme lists are correct
        self.time = 0
        # create a time-based seed and save it, and initialize our random generator with this seed
        self.seed = random.randrange(sys.maxsize)
        self.rng = np.random.default_rng(self.seed)
        # This option indicates if you want to include the specific sites for non-specific DNA binding proteins
        self.write_nonspecific_sites = False  # TODO: add this option as input
        # TODO: Create list of warnings in the future, and remove duplicates to print at the end of simulation
        # self.warnings = []  # List with warningss

        # -----------------------------------------------
        # We add new DNA sites which is the ones that we will link topos binding
        # Note: In the future there might be cases in which new enzymes will be added to the environment, so maybe,
        # these new sites will need to be created
        if topoisomerase_model == 'stochastic':
            topo_list = [environment for environment in self.environmental_list
                         if environment.enzyme_type == 'topo' or environment.enzyme_type == 'topoisomerase']
            for topo in topo_list:
                # The idea is that the topos will recognize these specific sites.
                # These sites will be created dynamically, and will have different names :0,1,2,3,...
                # I will have to ignore this first site in specific, because this one will be the overall, and is the
                # one that is output in the sites_df.csv
                t_site = Site(s_type='DNA_' + topo.name, name='DNA_' + topo.name + '_global',
                              start=1, end=self.size, k_min=0, k_max=0,
                              s_model=topo.binding_model + '_' + topo.name, oparams=topo.binding_oparams)
                self.site_list.append(t_site)

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

            # If we are modeling the topoisomerase binding as a stochastic process, we need to define the sites.
            # if self.topoisomerase_model == 'stochastic':
            #    self.define_topoisomerase_binding_sites()
            # BINDING
            # --------------------------------------------------------------
            # Apply binding model and get list of new enzymes
            new_enzyme_list = bm.binding_model(self.enzyme_list, self.environmental_list, self.dt, self.rng)
            # These new enzymes are lacking twist and superhelical, we need to fix them and actually add them
            # to the enzyme_list
            # But before, add the binding events to the log  (it's easier to do it first)
            #            self.add_binding_events_to_log(new_enzyme_list)
            self.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

            # EFFECT
            # --------------------------------------------------------------
            effects_list = em.effect_model(self.enzyme_list, self.environmental_list, self.dt,
                                           self.topoisomerase_model, self.mechanical_model)
            self.apply_effects(effects_list)

            # UNBINDING
            # --------------------------------------------------------------
            drop_list_index, drop_list_enzyme = bm.unbinding_model(self.enzyme_list, self.dt, self.rng)
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
            new_enzyme_list = bm.binding_model(self.enzyme_list, self.environmental_list, self.dt,
                                               self.rng)
            self.add_new_enzymes(new_enzyme_list)  # It also calculates fixes the twists and updates supercoiling

            # EFFECT
            # --------------------------------------------------------------
            effects_list = em.effect_model(self.enzyme_list, self.environmental_list, self.dt,
                                           self.topoisomerase_model, self.mechanical_model)
            self.apply_effects(effects_list)

            # UNBINDING
            # --------------------------------------------------------------
            drop_list_index, drop_list_enzyme = bm.unbinding_model(self.enzyme_list, self.dt,
                                                                   self.rng)
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
            d = {'type': enzyme.enzyme_type, 'name': enzyme.name, 'site': enzyme.site.name,
                 'position': enzyme.position, 'direction': enzyme.direction, 'size': enzyme.size, 'twist': enzyme.twist,
                 'superhelical': enzyme.superhelical}
            enzyme_aux.append(d)
        my_df = pd.DataFrame.from_dict(enzyme_aux)
        return my_df

    # Returns environmental list in the form of dataframe. This function is with the intention of outputting the system
    def environmental_list_to_df(self):
        environmental_aux = []  # This will be a list of dicts
        for environmental in self.environmental_list:
            d = {'type': environmental.enzyme_type, 'name': environmental.name,
                 'site_type': environmental.site_type, 'concentration': environmental.concentration,
                 'k_on': environmental.k_on, 'k_off': environmental.k_off, 'k_cat': environmental.k_cat,
                 'size': environmental.size}
            environmental_aux.append(d)
        my_df = pd.DataFrame.from_dict(environmental_aux)
        return my_df

    # Returns list of sites in the form of dataframe. This function is with the intention of outputting the system
    def site_list_to_df(self):
        site_aux = []  # This will be a list of dicts
        for site in self.site_list:
            d = {'type': site.site_type, 'name': site.name, 'start': site.start, 'end': site.end, 'k_min': site.k_min,
                 'k_max': site.k_max, 'model': site.site_model, 'oparams': site.oparams}
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
            # This is for enzymes that bind bare DNA
            if 'DNA' in site.site_type and '_global' in site.name:
                site_superhelical = self.superhelical
                site_twist = self.twist
            else:
                enzyme_before = [enzyme for enzyme in self.enzyme_list if enzyme.position <= site.start]
                if len(enzyme_before) <= 0:
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

            # This is for enzymes that bind bare DNA - Let's initialize it
            if 'DNA' in site.site_type and '_global' in site.name:
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
        self.site_list.append(Site(s_type='EXT', name='EXT', start=1, end=self.size, k_min=0, k_max=0, s_model=None,
                                   oparams=None))  # twist=self.twist, superhelical=self.superhelical) )

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
                    position_left, position_right = em.get_start_end_c(self.enzyme_list[0], self.enzyme_list[-1],
                                                                       self.size)

                else:  # For linear DNA
                    position_left = 0
                    position_right = self.size + 1

            else:  # If nothing is bound
                position_left = 0
                position_right = self.size + 1  # it is the same in this case for either linear or circular

            # TODO: I can distribute the supercoiling when defining the sites?
            #            for site in self.site_list:  # Distribute supercoiling -
            #                site.superhelical = self.superhelical
            for enzyme in self.enzyme_list:  # Distribute supercoiling -
                enzyme.superhelical = self.superhelical

            # Let's treat the boundaries of our system as objects.
            # ----------------------------------------------------------------
            # Notice that we don't specify the end
            # TODO: the site needs to be defined first
            extra_left = Enzyme(e_type='EXT', name='EXT_L', site=self.site_match('EXT', 'EXT'), position=position_left,
                                size=0, k_off=0, twist=0, superhelical=self.superhelical,
                                effect_model=None, effect_oparams=None)
            extra_right = Enzyme(e_type='EXT', name='EXT_R', site=self.site_match('EXT', 'EXT'),
                                 position=position_right,
                                 size=0, k_off=0, twist=0, superhelical=self.superhelical,
                                 effect_model=None, effect_oparams=None)

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
            enzyme.twist = em.calculate_twist(enzyme, self.enzyme_list[i + 1])

    # Updates the supercoiling/superhelical in enzymes
    def update_supercoiling(self):
        for i, enzyme in enumerate(self.enzyme_list[:-1]):
            enzyme.superhelical = em.calculate_supercoiling(enzyme, self.enzyme_list[i + 1])

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
        # (Maybe it shouldn't be provided)

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
            region_length = em.calculate_length(enzyme_before, enzyme_after)

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
                position_left, position_right = em.get_start_end_c(self.enzyme_list[1], self.enzyme_list[-2],
                                                                   self.size)
                self.enzyme_list[0].position = position_left
                self.enzyme_list[-1].position = position_right

            # We are still missing the supercoiling density and the excess of twist...
            # We need to partition the twist, so it is conserved...
            # --------------------------------------------------------
            # These quantities are the sizes of the new local domains on the left and right of the new enzyme
            new_length_left = em.calculate_length(enzyme_before, new_enzyme)
            new_length_right = em.calculate_length(new_enzyme, enzyme_after)

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
                    em.get_start_end_c(self.enzyme_list[1], self.enzyme_list[-2], self.size)
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
            position_left, position_right = em.get_start_end_c(self.enzyme_list[1],
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

            if enzyme.name == 'RNAP':
                size = abs(enzyme.site.start - enzyme.site.end + 1)
                output_environment = Environment(e_type='mRNA', name=enzyme.site.name, site_list=[], concentration=1,
                                                 k_on=0, k_off=0, k_cat=0, size=size, site_type=None, oparams=None)
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
                length = em.calculate_length(enzyme, next_enzyme)
                n_sites = int(length / topo.size)
                for n in range(n_sites):  # The 1+n is to leave some space 1 empty space between enzymes
                    start = enzyme.position + enzyme.size + topo.size * n + (1 + n)
                    end = 1 + enzyme.position + enzyme.size + topo.size * (1 + n)
                    if end > next_enzyme.position:  # Little break to avoid the overlapping of enzymes
                        continue
                    topo_site = Site(s_type='DNA_' + topo.name, name=str(s), start=start, end=end, k_min=0, k_max=0,
                                     s_model='stochastic_' + topo.name, oparams=None)
                    self.site_list.append(topo_site)

                    s = s + 1
                    topo.site_list.append(topo_site)
        self.sort_site_list()
        return

    # This function defines the binding sites of enzymes that recognize bare DNA, that means just DNA.
    # It partitions the DNA in N binding sites of size enzyme.size
    def define_bare_DNA_binding_sites(self):
        environment_list = [environment for environment in self.environmental_list
                            if environment.site_type == 'DNA']

        for environment in environment_list:
            # No point in defining topoisomerase binding sites if it's a continuum model
            if environment.enzyme_type == 'topo' and self.topoisomerase_model == 'continuum':
                continue
            n_sites = int(self.size / environment.size)
            s = 0
            for n in range(n_sites):
                start = 1 + environment.size * n
                end = environment.size * (1 + n)
                if end > self.size:  # Little break to avoid making it bigger than the actual plasmid
                    continue
                environment_site = Site(s_type='DNA_' + environment.name,
                                        name=str(s),
                                        start=start, end=end, k_min=0, k_max=0,
                                        s_model=environment.binding_model + '_' + environment.name,
                                        oparams=environment.binding_oparams)
                self.site_list.append(environment_site)

                s = s + 1
                environment.site_list.append(environment_site)

            # The next line makes the environmental recognize the specific binding site
            environment.site_type = 'DNA_' + environment.name

        self.sort_site_list()
        return

    def check_object_inputs(self):

        """
         Checks that the input objects (sites, environments, enzymes) have the correct parameters or if they were not
         provided, it assigns them to the default models if necessary
         """

        # TODO: Define list of available and supported models in params.
        #  If the model that is specified is not supported, then remove the model by model = None, and oparams=None
        # Sites
        # ================================
        for site in self.site_list:
            if site.site_model == 'maxmin' and site.oparams is None:
                site.site_model = 'sam'
                print('For site ', site.name, 'maxmin model selected but no parameters provided. '
                                              'Switching to sam model')

        # Environment
        # ================================
        for environment in self.environmental_list:

            # For the binding model
            # -------------------------------------------------------------------
            if environment.binding_oparams is None:

                # In case of topoisomerase enzyme and stochastic model
                if self.topoisomerase_model == 'stochastic' and 'topo' in environment.enzyme_type:
                    environment.binding_model = 'recognition'  # The default recognition model in case of stochastic
                    if 'topoI' in environment.name:
                        environment.binding_oparams = {'width': params.topo_b_w, 'threshold': params.topo_b_t}
                    elif 'gyrase' in environment.name:
                        environment.binding_oparams = {'width': params.gyra_b_w, 'threshold': params.gyra_b_t}

            # In case it is a continuum model, then there's no binding
            if self.topoisomerase_model == 'continuum' and 'topo' in environment.enzyme_type:
                environment.binding_model = None
                environment.binding_oparams = None

            # For the effect model
            # -------------------------------------------------------------------
            if environment.effect_oparams is None:
                # In case of topo
                if 'topo' in environment.enzyme_type:
                    if 'topoI' in environment.name:
                        environment.effect_oparams = {'k_cat': params.topo_k, 'width': params.topo_e_w,
                                                      'threshold': params.topo_e_t}
                    elif 'gyrase' in environment.name:
                        environment.effect_oparams = {'k_cat': params.gyra_k, 'width': params.gyra_e_w,
                                                      'threshold': params.gyra_e_t}
                # In case of RNAP
                elif 'RNAP' in environment.enzyme_type:
                    environment.effect_oparams = {'velocity': params.v0, 'gamma': params.gamma}

        # Enzyme
        # ================================
        for enzyme in self.enzyme_list:
            if enzyme.effect_oparams is None:
                # In case of topo
                if 'topo' in environment.enzyme_type:
                    if self.topoisomerase_model == 'continuum':
                        enzyme.effect_model = 'continuum'
                    if 'topoI' in environment.name:
                        environment.effect_oparams = {'k_cat': params.topo_k, 'width': params.topo_e_w,
                                                      'threshold': params.topo_e_t}
                    elif 'gyrase' in environment.name:
                        environment.effect_oparams = {'k_cat': params.gyra_k, 'width': params.gyra_e_w,
                                                      'threshold': params.gyra_e_t}
                # In case of RNAP
                elif 'RNAP' in environment.enzyme_type:
                    environment.effect_oparams = {'velocity': params.v0, 'gamma': params.gamma}
