import pandas as pd
from TORCphysics import binding_model as bm
import sys


class Site:
    """
    Site class.
    """

    def __init__(self, s_type, name, start, end, k_min, k_max, s_model_name, oparams_file, binding_model=None):
        """
        Initialize an instance of Site clase.

        Args:
            s_type (str): The syte type, e.g. gene. You can define overlapping sites, such as gene_TF, which would
            be a gene recognized by transcription factors TF. This particular site would have different functionality.
            name (str): The name of the site, e.g. tetA
            start (int): The starting position of the site.
            end (int): The ending position of the site.
            k_min (float): The minimum binding rate.
            k_max (float): The maximum binding rate.
            s_model_name (str): The name of the site model or binding model.
            oparams (str): Path to a dictionary of additional parameters relevant to the site_model s_model.

        Example:
            tetA_site = Site(
                s_type='gene',
                name='tetA',
                start=100,
                end=200,
                k_min=0.5,
                k_max=2.0,
                s_model='poisson',
                oparams=None
            )
        """
        self.oparams_file = oparams_file
        self.site_type = s_type
        self.name = name
        self.start = start
        self.end = end
        self.k_min = k_min
        self.k_max = k_max
        self.site_model = s_model_name
        self.direction = self.get_direction()
        self.
        self.oparams = oparams
        if isinstance(oparams, dict):
            self.oparams = oparams
            # print('is dict', self.name, self.site_type)
        else:
            # print('is not dict')
            self.oparams = None
            self.read_oparams(oparams)

    def read_oparams(self, oparams):
        """
        reads oparams that are related to the site_model
        """

        if oparams is None or oparams == 'none':
            self.oparams = None
        else:
            self.oparams = pd.read_csv(oparams).to_dict()

    def get_direction(self):
        """
        According start and end, gets the direction of the site.
        Only genes can have directions
        """

        direction = 0
        # Doesn't have direction if it's not a gene
        if self.site_type != 'gene':
            return direction
        else:
            if self.start < self.end:
                direction = 1
                return direction
            elif self.start > self.end:
                direction = -1
                return direction
            else:
                print("Cannot work out gene's direction")
                sys.exit()


class SiteFactory:
    """
    SiteFactory class. It is used to read csv input files and according to each entry, create a Site object
    """

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
                            k_min=float(row['k_min']), k_max=float(row['k_max']), s_model_name=row['model'],
                            oparams=row['oparams'])  # , twist=row['twist'], superhelical=['superhelical'])
            self.site_list.append(new_site)
