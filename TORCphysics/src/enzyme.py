import pandas as pd


class Enzyme:

    def __init__(self, e_type, name, site, position, size, k_cat, k_off, twist, superhelical):
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
        self.k_cat = k_cat
        self.k_off = k_off
        #I'm not sure if defining twist_front/behind was a good idea
        #self.twist_front = 0.0
        #self.twist_behind = 0.0
        #self.superhelical_front = 0.0
        #self.superhelical_behind = 0.0


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
                                position=float(row['position']), size=float(row['size']), k_cat=float(row['k_cat']),
                                k_off=float(row['k_off']),
                                twist=float(row['twist']), superhelical=float(row['superhelical']))

            self.enzyme_list.append(new_enzyme)

    def site_match(self, label):
        if label in [site.name for site in self.site_list]:
            # TODO check if this works!
            for site in self.site_list:
                if site.name == label:
                    return site  # the first one?
        else:
            return None
