import pandas as pd
import sys
from TORCphysics import params


class Environment:

    def __init__(self, e_type, name, site_list, concentration, k_on, k_off, k_cat, size, site_type,
                 oparams):  # , eff_model):
        self.oparams = None
        self.enzyme_type = e_type
        self.name = name
        self.site_list = site_list  # It recognizes a list of sites, rather than a specific site
        self.site_type = site_type  # We need to remember the type
        self.concentration = concentration
        self.k_on = k_on
        self.k_off = k_off
        self.k_cat = k_cat
        self.size = size
        self.read_oparams(oparams)

    def read_oparams(self, oparams):
        if oparams is None or oparams == 'none':
            # TODO: In the future would we need the topo model?
            #  Maybe if we want calibrated both continuum and stochastic
            if self.enzyme_type == 'topo' and ('topoI' in self.name):
                self.oparams = {'width': params.topo_w, 'threshold': params.topo_t}
            elif self.enzyme_type == 'topo' and ('gyrase' in self.name):
                self.oparams = {'width': params.gyra_w, 'threshold': params.gyra_t}
            else:
                self.oparams = None
        else:
            # TODO: Why it looks different when I load from the csv and causes errors when doing operations?
            #  Needs some testing! It certanly looks different than loading a dictionary. Maybe we need to make
            #  it a dictionary here
            self.oparams = pd.read_csv(oparams).to_dict()
            if self.enzyme_type == 'topo':
                if 'width' not in self.oparams or 'threshold' not in self.oparams:
                    print('Please provide width and threshold for topoisomerase parameters')
                    sys.exit()


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
                                          k_off=float(row['k_off']), k_cat=float(row['k_cat']),
                                          size=float(row['size']), site_type=row['site_type'],
                                          oparams=row['oparams'])  # eff_model=row['model'])
            self.environment_list.append(new_environment)

    def site_match(self, label):
        #        enzyme_before = [enzyme.position for enzyme in enzyme_list if enzyme.position <= site.start][-1]
        site_list = [site for site in self.site_list if site.site_type == label]
        return site_list

#        if label in [site.name for site in self.site_list]:
#            # TODO check if this works!
#            for site in self.site_list:
#                if site.name == label:
#                    return site  # the first one?
#        else:
#            return None
