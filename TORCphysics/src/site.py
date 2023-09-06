import pandas as pd
from TORCphysics import binding_model as bm
import sys


# TODO: It seems that the Site is working now! The next steps are as follows:
#  1.- Document Site and Environment.
#  2.- Fix Enzyme.
#  3.- Document Enzyme. 3.1- Test Enzyme. 3.2- Test Binding models. 3.3- Test Effect models
#  4.- Make the code works with the new changes, and don't forget to add the effective_size.
#  5.- Document binding/effect.
#  6.- Document everything

class Site:
    """    A class used to represent functioning sites (sequence) on the DNA.
    Attributes
    ----------
    site_type : str
     The syte type, e.g. gene. You can define overlapping sites, such as gene_TF, which would
     be a gene recognized by transcription factors TF. This particular site would have different functionality.
    name : str
     The name of the site, e.g. tetA
    start : float
     The starting position of the site.
    end : float
     The ending position of the site.
    k_on : float
     The minimum binding rate.
    binding_model_name : str, optional
     The name of the site model or binding model. It indicates how environmentals
     (enzymes on the environment) will bind the site
    binding_oparams_file : str, optional
     Path to a csv file with additional parameters relevant to the binding_model.
    binding_model : BindingModel class (from binding_model.py), optional
     The preloaded binding model to use. This binding model already contains any additional oparams parameters.
    binding_model_oparams : dict, optional
     Dictionary with parameters to include in the binding model.

    Notes
    ----------------
    Additional parameters (oparams) must be compatible with the given binding model, so do not create parameters that
    the model won't consider.
    When these oparams are not given, the code will load default parameters according the binding model used.
    """

    def __init__(self, site_type, name, start, end, k_on,
                 binding_model_name=None, binding_oparams_file=None,
                 binding_model=None, binding_model_oparams=None):
        """
        Initialize an instance of Site class.

        Args:
        site_type : str
         The syte type, e.g. gene. You can define overlapping sites, such as gene_TF, which would
         be a gene recognized by transcription factors TF. This particular site would have different functionality.
        name : str
         The name of the site, e.g. tetA
        start : float
         The starting position of the site.
        end : float
         The ending position of the site.
        k_on : float
         The minimum binding rate.
        binding_model_name : str, optional
         The name of the site model or binding model. It indicates how environmentals
         (enzymes on the environment) will bind the site
        binding_oparams_file : str, optional
         Path to a csv file with additional parameters relevant to the binding_model.
        binding_model : BindingModel class (from binding_model.py), optional
         The preloaded binding model to use. This binding model already contains any additional oparams parameters.
        binding_model_oparams : dict, optional
         Dictionary with parameters to include in the binding model.

        Example:
            tetA_site = Site(
                site_type='gene',
                name='tetA',
                start=100,
                end=200,
                k_min=0.5,
                k_max=2.0,
                binding_model_name='PoissonBinding'
            )
        """
        # Assign parameters
        self.site_type = site_type
        self.name = name
        self.start = start
        self.end = end
        self.k_on = k_on

        # Models
        self.binding_model_name = binding_model_name
        self.binding_oparams_file = binding_oparams_file
        self.binding_model = binding_model
        self.binding_model_oparams = binding_model_oparams

        self.check_inputs()

        self.direction = self.get_direction()

        self.get_models()  # Loads the binding model

    def check_inputs(self):

        if not isinstance(self.site_type, str) or self.site_type == '':
            raise ValueError('Error, sites must have a site type')
        if not isinstance(self.name, str) or self.name == '':
            raise ValueError('Error, sites must have a name')
        if not isinstance(self.start, float) and not isinstance(self.start, int):
            raise ValueError('Error, site start need a number')
        if not isinstance(self.end, float) and not isinstance(self.end, int):
            raise ValueError('Error, site end need a number')
        if not isinstance(self.k_on, float) and not isinstance(self.k_on, int):
            raise ValueError('Error, site k_on must be a number')
        #        if not isinstance(self.k_min, float) and not isinstance(self.k_min, int):
        #            raise ValueError('Error, site k_min must be a number')
        #        if not isinstance(self.k_max, float) and not isinstance(self.k_max, int):
        #            raise ValueError('Error, site k_max must be a number')

        # Binding model
        if self.binding_model_name == '' or self.binding_model_name == 'None' or self.binding_model_name == 'none':
            self.binding_model_name = None
        if self.binding_model == '' or self.binding_model == 'None' or self.binding_model == 'none':
            self.binding_model = None
        if (self.binding_oparams_file == '' or self.binding_oparams_file == 'None'
                or self.binding_oparams_file == 'none'):
            self.binding_oparams_file = None
        if (self.binding_model_oparams == '' or self.binding_model_oparams == 'None' or
                self.binding_model_oparams == 'none'):
            self.binding_model_oparams = None

        # Check k_on coincides with the one in oparams
        if self.binding_model_oparams is None:  # No dict
            if self.binding_oparams_file is None:  # No path
                self.binding_model_oparams = {'k_on', self.k_on}  # Then add k_on
        else:  # With dict
            if 'k_on' not in self.binding_model_oparams:  # But no k_on
                self.binding_model_oparams['k_on'] = self.k_on  # Add k_on - just in case
            else:  # There is k_on, but the site also has a k_on. Let's prioritise the one on the site.
                self.binding_model_oparams['k_on'] = self.k_on

    def get_models(self):

        # Binding Model
        self.binding_model, self.binding_model_name, self.binding_oparams_file, self.binding_model_oparams = (
            bm.get_binding_model(self.name, self.binding_model, self.binding_model_name,
                                 self.binding_oparams_file, self.binding_model_oparams))

        # Check one more that k_on is the same as the one in oparams.
        if self.binding_model is not None:
            self.binding_model_oparams['k_on'] = self.k_on  # Add k_on - just in case
            self.binding_model.k_on = self.k_on

    #            if 'k_on' not in self.binding_model_oparams:  # But no k_on
    #                self.binding_model_oparams['k_on'] = self.k_on  # Add k_on - just in case
    #                self.binding_model.k_on = self.k_on
    #            else:  # There is k_on, but the site also has a k_on. Let's prioritise the one on the site.
    #                self.binding_model_oparams['k_on'] = self.k_on
    #                self.binding_model.k_on = self.k_on

    #    def read_oparams(self, oparams):
    #        """
    #        reads oparams that are related to the site_model
    #        """##

    #        if oparams is None or oparams == 'none':
    #            self.oparams = None
    #        else:
    #            self.oparams = pd.read_csv(oparams).to_dict()

    def get_direction(self):
        """
        According start and end, gets the direction of the site.
        Only genes can have directions
        """

        direction = 0
        # Doesn't have direction if it's not a gene
        if self.site_type.lower() != 'gene':
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

    def __init__(self, filename=None):
        self.filename = filename
        self.site_list = []
        if filename:
            self.read_csv()

    def get_site_list(self):
        return self.site_list

    def read_csv(self):
        df = pd.read_csv(self.filename)
        for index, row in df.iterrows():
            new_site = Site(site_type=row['type'], name=row['name'], start=float(row['start']), end=float(row['end']),
                            k_on=float(row['k_on']), binding_model_name=row['binding_model'],
                            binding_oparams_file=row['binding_oparams'])
            # new_site = Site(site_type=row['type'], name=row['name'], start=float(row['start']), end=float(row['end']),
            #  k_min=float(row['k_min']), k_max=float(row['k_max']), binding_model_name=row['model'],
            #                            binding_model_oparams=row['oparams'])
            self.site_list.append(new_site)
