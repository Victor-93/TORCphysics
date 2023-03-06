import pandas as pd
import random
import sys
from TORCphysics import Site, SiteFactory, Enzyme, EnzymeFactory, EnvironmentFactory#, mechanical_model
from TORCphysics import mechanical_model as mm


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
        self.size = 0
        self.twist = None
        self.superhelical = None
        self.read_csv()  # Here, it gets the name,structure, etc
        self.site_list = SiteFactory(sites_filename).site_list
        self.enzyme_list = EnzymeFactory(enzymes_filename, self.site_list).enzyme_list
        self.environmental_list = EnvironmentFactory(environment_filename, self.site_list).environment_list
        self.time = 0
        # create a time-based seed and save it
        self.seed = random.randrange(sys.maxsize)

        #-----------------------------------------------
        #We add another SIDD which is the one we will link topos binding with DNA
        self.site_list.append( Site(s_type='DNA', name='DNA', start=1, end=self.size, k_min=0, k_max=0, s_model=None,
                                    oparams=None))
        #Sort list of enzymes and sites by position/start
        self.sort_lists()
        # Distribute twist/supercoiling
        self.add_fake_boundaries()
        self.sort_lists()
        # Define local sites
#        self.calculate_local_sites()

        # TODO: I have to think about the update, what does the update do? Do I really need one? Or maybe many updates
        #  for example, a twist update, a supercoiling update, a local sites update?
        #  maybe I can update sites as well,
        # TODO: My update should update twist, supercoiling and the local sites?

        # TODO: sites, enzymes and superhelical have been loaded, but still needs to define the local sites/domains
        #  (bare DNA) - or is it the update?

    def sort_lists(self):
        self.enzyme_list.sort(key=lambda x: x.position)
        self.site_list.sort(key=lambda x: x.start)
#    def calculate_local_sites(self):
#    #This function calculates the local sites with naked DNA.
#    #If there are N enzymes, then there is N-1 local sites (or local domains).
#        for i in range(self.get_num_enzymes()):
#            self.site_list.append( Site())
#        for enzyme in self.enzyme_list:
#            print(0)

    def add_fake_boundaries(self):
        #I need to add a fake site so I can link the fake boundaries
        self.site_list.append( Site(s_type='EXT', name='EXT', start=1, end=self.size, k_min=0, k_max=0, s_model=None,
                                    oparams=None))#, twist=self.twist, superhelical=self.superhelical) )

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
                    position_left, position_right = mm.get_start_end_c(self.enzyme_list[0], self.enzyme_list[-1], self.size)

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
                #enzyme.superhelical_front = self.superhelical
                #enzyme.superhelical_behind = self.superhelical
                enzyme.superhelical = self.superhelical

            # Let's treat the boundaries of our system as objects.
            # ----------------------------------------------------------------
            # Notice that we don't specify the end
            # TODO: the site needs to be defined first
            extra_left = Enzyme(e_type='EXT', name='EXT_L', site=self.site_match('EXT'), position=position_left,
                                size=0, twist=0, superhelical=self.superhelical)
            extra_right = Enzyme(e_type='EXT', name='EXT_R', site=self.site_match('EXT'), position=position_right,
                                 size=0, twist=0, superhelical=self.superhelical)

            self.enzyme_list.append( extra_left)
            self.enzyme_list.append( extra_right)
            self.sort_lists()
            #And finally, update the twist
            self.update_twist()

            print(0)
            # WARNING!!!!
            # There could be a big mistake in case of linear structures that have a NAP in positions 1 or nbp

    #This functions updates twist in enzymes
    def update_twist(self):

        for i, enzyme in enumerate(self.enzyme_list[:-1]):
            enzyme.twist = mm.calculate_twist(enzyme, self.enzyme_list[i+1])

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
