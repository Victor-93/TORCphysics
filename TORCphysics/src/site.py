import pandas

class site:

    def __init__(self, s_type, name, start, end, k_min, k_max, s_model, oparams):
        self.site_type = s_type
        self.name = name
        self.start = start
        self.end = end
        self.k_min = k_min
        self.k_max = k_max
        self.site_model = s_model
        self.oparams = self.set_up_oparams(oparams)

    def set_up_oparams(self, param_string):
        oparams = {}
        if self.site_type == x:
            # TODO: process param_string for different site types
            oparams[p] = x.split(",")[3]
        return oparams


class SiteFactory:

    def __init__(self, filename):
        self.filename = filename
        self.site_list = []
        self.read_csv()

    def get_site_list(self):
        return self.site_list

    def read_csv(self):
        df = pandas.read_csv()
