import pandas as pd
import random
import sys
from TORCphysics import Site, SiteFactory, Enzyme, EnzymeFactory, EnvironmentFactory, Event, Log, params
from TORCphysics import effect_model as em
from TORCphysics import binding_model as bm


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
        self.time = 0
        # create a time-based seed and save it
        self.seed = random.randrange(sys.maxsize)

        # -----------------------------------------------
        # We add another SIDD which is the one we will link topos binding with DNA
        self.site_list.append(Site(s_type='DNA', name='DNA', start=1, end=self.size, k_min=0, k_max=0, s_model=None,
                                   oparams=None))
        # Sort list of enzymes and sites by position/start
        self.sort_lists()
        # Distribute twist/supercoiling
        self.add_fake_boundaries()
        self.sort_lists()

        # TODO: test that the update of global superhelical correctly  matches the initial sigma. I think the problem is
        #  due to the sizes of enzymes. I have two options: 1-include those sizes, 2-include the twist
        #  of the bound region. For now I'll take option 2, when an enzyme binds, the region becomes flat(relax)
        self.update_global_twist()
        self.update_global_superhelical()

        # Let's initialize the log
        self.log = Log(self.size, self.frames, self.frames*self.dt, self.structure, self.name, self.site_list,
                       self.twist, self.superhelical)

    # Define local sites

    #        self.calculate_local_sites()

    # TODO: I have to think about the update, what does the update do? Do I really need one? Or maybe many updates
    #  for example, a twist update, a supercoiling update, a local sites update?
    #  maybe I can update sites as well,
    # TODO: My update should update twist, supercoiling and the local sites?

    # TODO: sites, enzymes and superhelical have been loaded, but still needs to define the local sites/domains
    #  (bare DNA) - or is it the update?

    # This one runs the simulation
    # TODO: finish the run
    # TODO: test all these functions
    # TODO: Think about the updating of twist/supercoiling. When new enzymes bind, they can have an impact on the
    #  supercoiling, in which depending on the size of the local domain, the supercoiling can increase considerably.
    #  So when, these enzymes are added, they start moving?  Or is it assumed that during this dt it took them a lot
    #  of time to bind? And when the effects are taking place, the enzymes react to these new supercoiling/twist or
    #  they feel the one prior binding? Maybe if dt is sufficiently small, it just doesn't matter?
    def run(self):

        for frame in range(self.frames):
#            print('frame =', frame)
            self.frame = frame
            self.time = frame * self.dt
            # BINDING
            # --------------------------------------------------------------
            # Apply binding model and get list of new enzymes
            new_enzyme_list = bm.binding_model(self.enzyme_list, self.environmental_list, self.dt)
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
            drop_list = bm.unbinding_model(self.enzyme_list)
            self.drop_enzymes(drop_list)
#            self.add_unbinding_events_to_log(drop_list)

            # UPDATE GLOBALS
            # --------------------------------------------------------------
            self.update_global_twist()
            self.update_global_superhelical()

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
        self.superhelical = self.twist / (params.w0 * (self.size - sum(enzyme.size for enzyme in self.enzyme_list)))

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
                    position_left = 1
                    position_right = self.size

            else:  # If nothing is bound
                position_left = 1
                position_right = self.size  # it is the same in this case for either linear or circular

            # TODO: I can distribute the supercoiling when defining the sites?
            #            for site in self.site_list:  # Distribute supercoiling -
            #                site.superhelical = self.superhelical
            for enzyme in self.enzyme_list:  # Distribute supercoiling -
                # enzyme.superhelical_front = self.superhelical
                # enzyme.superhelical_behind = self.superhelical
                enzyme.superhelical = self.superhelical

            # Let's treat the boundaries of our system as objects.
            # ----------------------------------------------------------------
            # Notice that we don't specify the end
            # TODO: the site needs to be defined first
            extra_left = Enzyme(e_type='EXT', name='EXT_L', site=self.site_match('EXT'), position=position_left,
                                size=0, twist=0, superhelical=self.superhelical)
            extra_right = Enzyme(e_type='EXT', name='EXT_R', site=self.site_match('EXT'), position=position_right,
                                 size=0, twist=0, superhelical=self.superhelical)

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

    def site_match(self, label):
        if label in [site.name for site in self.site_list]:
            # TODO check if this works!
            for site in self.site_list:
                if site.name == label:
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

            # Before updating local parameters, create the new event and add it to log
            # --------------------------------------------------------------------------
            new_event = Event(self.time, self.frame, 'binding_event', enzyme_before.twist,
                              enzyme_before.superhelical,
                              self.twist, self.superhelical, new_enzyme.site, new_enzyme, new_enzyme.position)
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
            # is assumed to be relaxed by the enzyme
            new_twist_left = region_twist * ((new_length_left + 0.5 * new_enzyme.size) / region_length)
            new_twist_right = region_twist * ((new_length_right + 0.5 * new_enzyme.size) / region_length)

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
                    enzyme_before.twist = new_twist_left
                    self.enzyme_list[-1].twist = new_twist_left
                    new_enzyme.twist = new_twist_right

                # There is one EXT at the right
                elif enzyme_before.name != 'EXT_L' and enzyme_before.name == 'EXT_R':
                    enzyme_before.twist = new_twist_left
                    self.enzyme_list[0] = new_twist_right
                    self.enzyme_list[-1] = new_twist_right
                    new_enzyme.twist = new_twist_right

                # In any other case where there's no neighbour boundaries
                else:
                    enzyme_before.twist = new_twist_left
                    new_enzyme.twist = new_twist_right

            # ------------LINEAR DNA--------------------
            else:
                enzyme_before.twist = new_twist_left
                new_enzyme.twist = new_twist_right

            # Now add the enzyme to the list, sort it
            self.enzyme_list.append(new_enzyme)
            self.sort_enzyme_list()

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
                if self.enzyme_list[i - 1].name == 'EXT_L' and self.enzyme_list[i + 1].name != 'EXT_R':
                    self.enzyme_list[self.get_num_enzymes() - 2].twist += self.enzyme_list[i].twist
                    self.enzyme_list[0].twist = self.enzyme_list[self.get_num_enzymes() - 2].twist

                # ------.......E.......O.....---------- / Twist in O is added to E
                else:
                    self.enzyme_list[i - 1].twist += self.enzyme_list[i].twist

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
        if self.circle:
            if self.get_num_enzymes() > 2:
                self.enzyme_list[0].position, self.enzyme_list[-1].position = \
                    em.get_start_end_c(self.enzyme_list[1], self.enzyme_list[-2], self.size)
            else:
                self.enzyme_list[0].position = 1
                self.enzyme_list[-1].position = self.size

        self.sort_enzyme_list()
        self.update_supercoiling()

        # Now that the global supercoiling is updated, let's add the new unbinding events to the log
        for new_event in new_events:
            new_event.global_superhelical = self.superhelical
            new_event.global_twist = self.twist

            # And add it to the log
            self.log.metadata.append(new_event)

    # Apply effects in effects_list
    def apply_effects(self, effects_list):
        # And apply the effects for the specified enzymes in the effects_list
        for effect in effects_list:
            self.enzyme_list[effect.index].position += effect.position
            self.enzyme_list[effect.index].twist += effect.twist_right
            # In case we affect the boundary on the left - it affects the last (not fake) enzyme
            # because the fake boundaries don't exist and just reflect the first and last enzymes.
            if self.circle and effect.index == 1:
                self.enzyme_list[self.get_num_enzymes() - 2].twist += effect.twist_left
            else:  # In any other case just update the enzyme on the left
                self.enzyme_list[effect.index - 1].twist += effect.twist_left

        # And update supercoiling - because twist was modified
        self.update_supercoiling()

