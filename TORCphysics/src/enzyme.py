import pandas as pd
from TORCphysics import Site
from TORCphysics import effect_model as em
from TORCphysics import unbinding_model as ubm


class Enzyme:

    def __init__(self, e_type, name, site, position, size, effective_size, twist, superhelical,
                 effect_model_name=None, effect_oparams_file=None, effect_model=None, effect_model_oparams=None,
                 unbinding_model_name=None, unbinding_oparams_file=None, unbinding_model=None,
                 unbinding_model_oparams=None):

        # Assign parameters
        self.enzyme_type = e_type
        self.name = name
        self.site = site
        self.position = position
        self.size = size
        self.effective_size = effective_size
        self.twist = twist
        self.superhelical = superhelical

        # Assign models
        self.effect_model_name = effect_model_name
        self.effect_oparams_file = effect_oparams_file
        self.effect_model = effect_model
        self.effect_model_oparams = effect_model_oparams

        self.unbinding_model_name = unbinding_model_name
        self.unbinding_oparams_file = unbinding_oparams_file
        self.unbinding_model = unbinding_model
        self.unbinding_model_oparams = unbinding_model_oparams

        # Verify inputs
        self.check_inputs()

        # If it passed the inputs, then calculate the values associated with the site
        self.start = self.site.start
        self.end = self.site.end
        self.direction = self.site.get_direction()

        # Loads the binding, effect and unbinding models if given.
        self.get_models()

    def check_inputs(self):
        """ Checks that Enzyme parameters are of the correct type.
        """

        if not isinstance(self.enzyme_type, str) or self.enzyme_type == '':
            raise ValueError('Error, enzymes must have a type')
        if not isinstance(self.name, str) or self.name == '':
            raise ValueError('Error, enzymes must have a name')
        if not isinstance(self.site, Site):
            print('Error for enzyme ', self.name)
            raise ValueError('Error, (bound) Enzymes must be linked to a Site')
        if not isinstance(self.position, float) and not isinstance(self.position, int):
            raise ValueError('Error, enzymes need a number for their position')
        if not isinstance(self.size, float) and not isinstance(self.size, int):
            raise ValueError('Error, enzymes need a number for size')
        if not isinstance(self.effective_size, float) and not isinstance(self.effective_size, int):
            raise ValueError('Error, enzymes need a number for effective size')
#        if not isinstance(self.start, float) and not isinstance(self.start, int):
#            raise ValueError('Error, enzymes need a number for start')
#        if not isinstance(self.end, float) and not isinstance(self.end, int):
#            raise ValueError('Error, enzymes need a number for start')
#        if not isinstance(self.direction, float) and not isinstance(self.direction, int):
#            raise ValueError('Error, enzymes need a number for direction')
        if not isinstance(self.twist, float) and not isinstance(self.twist, int):
            raise ValueError('Error, enzymes need a number for twist')
        if not isinstance(self.superhelical, float) and not isinstance(self.superhelical, int):
            raise ValueError('Error, enzymes need a number for superhelical')

        if self.effective_size > self.size:
            print('Error, effective size effective_size cannot be larger than size')
            print('effective_size:', self.effective_size)
            print('size', self.size)
            print('For enzyme ', self.name)
            raise ValueError('Error: effective_size > size')

        # Effect model
        if self.effect_model_name == '' or self.effect_model_name == 'None' or self.effect_model_name == 'none':
            self.effect_model_name = None
        if self.effect_model == '' or self.effect_model == 'None' or self.effect_model == 'none':
            self.effect_model = None
        if self.effect_oparams_file == '' or self.effect_oparams_file == 'None' or self.effect_oparams_file == 'none':
            self.effect_oparams_file = None
        if (self.effect_model_oparams == '' or self.effect_model_oparams == 'None' or
                self.effect_model_oparams == 'none'):
            self.effect_model_oparams = None

        # Unbinding model
        if (self.unbinding_model_name == '' or self.unbinding_model_name == 'None'
                or self.unbinding_model_name == 'none'):
            self.unbinding_model_name = None
        if self.unbinding_model == '' or self.unbinding_model == 'None' or self.unbinding_model == 'none':
            self.unbinding_model = None
        if (self.unbinding_oparams_file == '' or self.unbinding_oparams_file == 'None'
                or self.unbinding_oparams_file == 'none'):
            self.unbinding_oparams_file = None
        if (self.unbinding_model_oparams == '' or self.unbinding_model_oparams == 'None' or
                self.unbinding_model_oparams == 'none'):
            self.unbinding_model_oparams = None

    def get_models(self):
        """ Loads the Enzyme's effect and unbinding models (if given).
        """

        # Effect Model
        self.effect_model, self.effect_model_name, self.effect_oparams_file, self.effect_model_oparams = (
            em.get_effect_model(self.name, self.effect_model, self.effect_model_name,
                                self.effect_oparams_file, self.effect_model_oparams))

        # Unbinding Model
        self.unbinding_model, self.unbinding_model_name, self.unbinding_oparams_file, self.unbinding_model_oparams = (
            ubm.get_unbinding_model(self.name, self.unbinding_model, self.unbinding_model_name,
                                    self.unbinding_oparams_file, self.unbinding_model_oparams))




class EnzymeFactory:

    def __init__(self, filename=None, site_list=None):
        self.filename = filename
        if site_list is not None:  # In case site_list is given but is not a list
            if not isinstance(site_list, list):
                raise ValueError('Error in EnzymeFactory. site_list must be a list if given.')
        self.site_list = site_list
        self.enzyme_list = []
        if filename:
            if site_list is None:  # In case site_list is given but is not a list
                raise ValueError('Error in EnzymeFactory. filename provided but site_list is missing.')
            if len(site_list) == 0:
                raise ValueError('Error in EnzymeFactory. filename provided but empty site_list.')
            self.read_csv()

    def get_enzyme_list(self):
        return self.enzyme_list

# In these inputs, site is actually giving the site's name. So there cannot be multiple names?
    def read_csv(self):
        df = pd.read_csv(self.filename)
        for index, row in df.iterrows():
            new_enzyme = Enzyme(e_type=row['type'], name=row['name'],
                                site=self.site_match(row['site']), position=float(row['position']),
                                size=float(row['size']), effective_size=float(row['effective_size']),
                                twist=float(row['twist']), superhelical=float(row['superhelical']),
                                effect_model_name=row['effect_model'],
                                effect_oparams_file=row['effect_oparams'],
                                unbinding_model_name=row['unbinding_model'],
                                unbinding_oparams_file=row['unbinding_oparams'])

            self.enzyme_list.append(new_enzyme)

    def site_match(self, label):
        if label in [site.name for site in self.site_list]:
            # TODO check if this works!
            for site in self.site_list:
                if site.name == label:
                    return site  # the first one?
        else:
            return None
