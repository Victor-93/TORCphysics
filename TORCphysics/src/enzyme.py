import pandas as pd
from TORCphysics import Site


class Enzyme:

    # TODO: Falta anadir el unbinding.
    def __init__(self, e_type, name, site, position, size, effective_size, k_off, twist, superhelical,
                 effect_model_name=None, effect_oparams_file=None, effect_model=None, effect_model_oparams=None):

        # Assign parameters
        self.enzyme_type = e_type
        self.name = name
        self.site = site
        self.position = position
        self.size = size
        self.effective_size = effective_size
        self.start = self.site.start
        self.end = self.site.end
        self.direction = self.site.get_direction()
        self.twist = twist
        self.superhelical = superhelical
        self.k_off = k_off

        # Assign models
        self.effect_model_name = effect_model_name
        self.effect_oparams_file = effect_oparams_file
        self.effect_model = effect_model
        self.effect_model_oparams = effect_model_oparams

        # Verify inputs
        self.check_inputs()

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
            raise ValueError('Error, (bound) Enzymes must be linked to a Site')
        if not isinstance(self.position, float) and not isinstance(self.position, int):
            raise ValueError('Error, enzymes need a number for their position')
        if not isinstance(self.size, float) and not isinstance(self.size, int):
            raise ValueError('Error, enzymes need a number for size')
        if not isinstance(self.effective_size, float) and not isinstance(self.effective_size, int):
            raise ValueError('Error, enzymes need a number for effective size')
        if not isinstance(self.start, float) and not isinstance(self.start, int):
            raise ValueError('Error, enzymes need a number for start')
        if not isinstance(self.end, float) and not isinstance(self.end, int):
            raise ValueError('Error, enzymes need a number for start')
        if not isinstance(self.direction, float) and not isinstance(self.direction, int):
            raise ValueError('Error, enzymes need a number for direction')
        if not isinstance(self.twist, float) and not isinstance(self.twist, int):
            raise ValueError('Error, enzymes need a number for twist')
        if not isinstance(self.superhelical, float) and not isinstance(self.superhelical, int):
            raise ValueError('Error, enzymes need a number for superhelical')
        if not isinstance(self.k_off, float) and not isinstance(self.k_off, int):
            raise ValueError('Error, enzymes need a number for k_off')

        if self.effective_size > self.size:
            print('Error, effective size effective_size cannot be larger than size')
            print('effective_size:', self.effective_size)
            print('size', self.size)
            print('For environmental ', self.name)
            raise ValueError('Error: effective_size > size')


    def read_oparams(self, oparams):
        """
        Reads oparams that are related to the site_model
        """

        if oparams is None or oparams == 'none':
            self.effect_oparams = None
        else:
            self.effect_oparams = pd.read_csv(oparams).to_dict()


class EnzymeFactory:

    def __init__(self, filename, site_list):
        self.filename = filename
        self.enzyme_list = []
        self.site_list = site_list
        self.read_csv()

    def get_enzyme_list(self):
        return self.enzyme_list

# In these inputs, site is actually giving the site's name. So there cannot be multiple names?
    def read_csv(self):
        df = pd.read_csv(self.filename)
        for index, row in df.iterrows():
            new_enzyme = Enzyme(e_type=row['type'], name=row['name'], site=self.site_match(row['site']),
                                position=float(row['position']), size=float(row['size']), k_off=float(row['k_off']),
                                twist=float(row['twist']), superhelical=float(row['superhelical']),
                                effect_model=row['effect_model'], effect_oparams=row['effect_oparams'])

            self.enzyme_list.append(new_enzyme)

    def site_match(self, label):
        if label in [site.name for site in self.site_list]:
            # TODO check if this works!
            for site in self.site_list:
                if site.name == label:
                    return site  # the first one?
        else:
            return None
