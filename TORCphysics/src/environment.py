import pandas as pd
import sys
from TORCphysics import params


class Environment:

    def __init__(self, e_type, name, site_list, concentration, k_on, k_off, size, site_type,
                 binding_model, binding_oparams, effect_model, effect_oparams):  # , eff_model):
        self.binding_oparams = None
        self.effect_oparams = None
        self.enzyme_type = e_type
        self.name = name
        self.site_list = site_list  # It recognizes a list of sites, rather than a specific site
        self.site_type = site_type  # We need to remember the type
        self.concentration = concentration
        self.k_on = k_on
        self.k_off = k_off
        self.size = size
        self.binding_model = binding_model
        self.effect_model = effect_model
        self.read_oparams(binding_oparams, effect_oparams)
        self.binding_oparams_file = binding_oparams
        self.effect_oparams_file = effect_oparams

    def read_oparams(self, binding_oparams, effect_oparams):

        # For the binding model
        # -------------------------------------------------------------------
        if binding_oparams is None or binding_oparams == 'none':
            self.binding_oparams = None
        else:
            self.binding_oparams = pd.read_csv(binding_oparams).to_dict()

        # For the effect model
        # -------------------------------------------------------------------
        if effect_oparams is None or effect_oparams == 'none':
            self.effect_oparams = None
        else:
            self.effect_oparams = pd.read_csv(effect_oparams).to_dict()




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
                                          k_off=float(row['k_off']), size=float(row['size']),
                                          site_type=row['site_type'],
                                          binding_model=row['binding_model'],
                                          binding_oparams=row['binding_oparams'],
                                          effect_model=row['effect_model'],
                                          effect_oparams=row['effect_oparams'])
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
