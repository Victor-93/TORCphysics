import pandas as pd


class Environment:

    def __init__(self, e_type, name, site_list, concentration, k_on, k_off, size, eff_model):
        self.enzyme_type = e_type
        self.name = name
        self.site_list = site_list  # It recognizes a list of sites, rather than a specific site
        self.concentration = concentration
        self.k_on = k_on
        self.k_off = k_off
        self.size = size
        self.eff_model = eff_model


class EnvironmentFactory:
    def __init__(self, filename, site_list):
        self.filename = filename
        self.environment_list = []
        self.site_list = site_list
        self.read_csv()

    def get_environment_list(self):
        return self.environment_list

    def read_csv(self):
        df = pd.read_csv(self.filename)
        for index, row in df.iterrows():
            new_environment = Environment(e_type=row['type'], name=row['name'],
                                          site_list=self.site_match(row['site_type']),
                                          concentration=float(row['concentration']), k_on=float(row['k_on']),
                                          k_off=float(row['k_off']), size=float(row['size']), eff_model=row['model'])
            self.environment_list.append(new_environment)

    def site_match(self, label):
        enzyme_before = [enzyme.position for enzyme in enzyme_list if enzyme.position <= site.start][-1]
        site_list = [site for site in self.site_list if site.site_type == label]
        return site_list

#        if label in [site.name for site in self.site_list]:
#            # TODO check if this works!
#            for site in self.site_list:
#                if site.name == label:
#                    return site  # the first one?
#        else:
#            return None
