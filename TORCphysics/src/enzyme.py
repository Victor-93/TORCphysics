import pandas as pd


class Enzyme:

    def __init__(self, e_type, name, site, position, size, k_off, twist, superhelical, effect_model, effect_oparams):
        self.effect_oparams = None
        self.enzyme_type = e_type
        self.name = name
        self.site = site
        self.position = position
        self.size = size
        self.start = self.site.start
        self.end = self.site.end
        self.direction = self.site.get_direction()
        self.twist = twist
        self.superhelical = superhelical
        self.k_off = k_off
        self.effect_model = effect_model
        self.read_oparams(effect_oparams)
        self.effect_oparams_file = effect_oparams

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
