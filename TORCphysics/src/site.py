import pandas as pd
import sys


class Site:

    def __init__(self, s_type, name, start, end, k_min, k_max, s_model, oparams):
        self.site_type = s_type
        self.name = name
        self.start = start
        self.end = end
        self.k_min = k_min
        self.k_max = k_max
        self.site_model = s_model
        self.oparams = oparams
        self.direction = self.get_direction()
#        self.oparams = self.set_up_oparams(oparams)
        # TODO: add direction

    def set_up_oparams(self, param_string):
        oparams = {}
        if self.site_type == x:
            # TODO: process param_string for different site types
            oparams[p] = x.split(",")[3]
        return oparams

    # According start and end, gets the direction.
    # Only genes can have directions
    def get_direction(self):
        direction = 0
        # Doesn't have direction if it's not a gene
        if self.site_type != 'gene':
            return direction
        else:
            if self.start < self.end:
                direction = -1
                return direction
            elif self.start > self.end:
                direction = 1
                return direction
            else:
                print("Cannot work out gene's direction")
                sys.exit()


class SiteFactory:

    def __init__(self, filename):
        self.filename = filename
        self.site_list = []
        self.read_csv()

    def get_site_list(self):
        return self.site_list

    def read_csv(self):
        df = pd.read_csv(self.filename)
        for index, row in df.iterrows():
            new_site = Site(s_type=row['type'], name=row['name'], start=float(row['start']), end=float(row['end']),
                            k_min=float(row['k_min']), k_max=float(row['k_max']), s_model=row['model'],
                            oparams=row['oparams'])
            self.site_list.append(new_site)
