import pandas as pd


class Enzyme:

    def __init__(self, e_type, name, site, position, size, twist, superhelical):
        self.enzyme_type = e_type
        self.name = name
        self.site = site
        self.position = position
        self.size = size
        self.twist = twist
        self.superhelical = superhelical
        self.start = self.site.start
        self.end = self.site.end
        self.direction = self.site.get_direction()


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
                                position=float(row['position']), size=float(row['size']), twist=float(row['twist']),
                                superhelical=float(row['superhelical']))

            self.enzyme_list.append(new_enzyme)

    def site_match(self, label):
        if label in [site.name for site in self.site_list]:
            # TODO check if this works!
            for site in self.site_list:
                if site.name == label:
                    return site  # the first one?
        else:
            return None
