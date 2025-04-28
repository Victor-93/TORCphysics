import pandas as pd
from TORCphysics import binding_model as bm


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
    global_site : bool
        Indicates if the site is a global site. Global sites are not actually bound by enzymes, but they keep track of
        the number of times an enzyme/molecule that recognizes bare DNA, bound the DNA anywhere.
    direction : int
        The direction of sites: -1 = left; 0 = no direction; 1 = right
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
    Additional parameters (oparams) must be compatible with the given binding, so please,
    do not create parameters that the models won't consider.
    When these oparams are not given, the code will load default parameters according the model indicated.
    These default parameters are saved in params.py
    """

    def __init__(self, site_type, name, start, end, k_on,
                 binding_model_name=None, binding_oparams_file=None,
                 binding_model=None, binding_model_oparams=None,
                 global_site=False):
        """ The constructor of  Site class.

        Parameters
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
        global_site : bool
            Indicates if the site is a global site. Global sites are not actually bound by enzymes, but they keep track of
            the number of times an enzyme/molecule that recognizes bare DNA, bound the DNA anywhere.


        Example
        ----------
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
        self.global_site = global_site

        # Assign models
        self.binding_model_name = binding_model_name
        self.binding_oparams_file = binding_oparams_file
        self.binding_model = binding_model
        self.binding_model_oparams = binding_model_oparams

        # Verify inputs
        self.check_inputs()

        # Get the direction
        self.direction = self.get_direction()

        # Loads the binding model
        self.get_models()

    def check_inputs(self):
        """ Checks that Site parameters are of the correct type.
        """

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
        if (self.binding_model_name == '' or self.binding_model_name == 'None' or self.binding_model_name == 'none'
                or self.binding_model_name == 'nan'):
            self.binding_model_name = None
        if (self.binding_model == '' or self.binding_model == 'None' or self.binding_model == 'none' or
                self.binding_model == 'nan'):
            self.binding_model = None
        if (self.binding_oparams_file == '' or self.binding_oparams_file == 'None'
                or self.binding_oparams_file == 'none' or self.binding_oparams_file == 'nan'):
            self.binding_oparams_file = None
        if (self.binding_model_oparams == '' or self.binding_model_oparams == 'None' or
                self.binding_model_oparams == 'none' or self.binding_model_oparams == 'nan'):
            self.binding_model_oparams = None

        # Check k_on coincides with the one in oparams
        if self.binding_model_oparams is None:  # No dict
            if self.binding_oparams_file is None:  # No path
                self.binding_model_oparams = {'k_on', self.k_on}  # Then add k_on
        else:  # With dict
            self.binding_model_oparams['k_on'] = self.k_on  # Add k_on - just in case
            #  if 'k_on' not in self.binding_model_oparams:  # But no k_on
            #    self.binding_model_oparams['k_on'] = self.k_on  # Add k_on - just in case
            #  else:  # There is k_on, but the site also has a k_on. Let's prioritise the one on the site.
            #    self.binding_model_oparams['k_on'] = self.k_on

    def get_models(self):
        """ Loads the Site's binding model (if given).
        """

        # Binding Model
        # print(self.name)
        self.binding_model, self.binding_model_name, self.binding_oparams_file, self.binding_model_oparams = (
            bm.get_binding_model(self.name, self.binding_model, self.binding_model_name,
                                 self.binding_oparams_file, self.binding_model_oparams))

        # Check one more that k_on is the same as the one in oparams.
        if self.binding_model is not None:
            self.binding_model_oparams['k_on'] = self.k_on  # Add k_on - just in case
            self.binding_model.k_on = self.k_on
            if hasattr(self.binding_model, 'k_max'):  # This is for some models that have maximum rate
                self.binding_model_oparams['k_max'] = self.k_on  # Add k_on - just in case
                self.binding_model.k_max = self.k_on

    #            if 'k_on' not in self.binding_model_oparams:  # But no k_on
    #                self.binding_model_oparams['k_on'] = self.k_on  # Add k_on - just in case
    #                self.binding_model.k_on = self.k_on
    #            else:  # There is k_on, but the site also has a k_on. Let's prioritise the one on the site.
    #                self.binding_model_oparams['k_on'] = self.k_on
    #                self.binding_model.k_on = self.k_on

    def get_direction(self):
        """
        According start and end, gets the direction of the site.
        Only genes can have directions
        """

        direction = 0
        # Doesn't have direction if it's not a gene
        if 'gene' in self.site_type.lower():
            if self.start < self.end:
                direction = 1
                return direction
            elif self.start > self.end:
                direction = -1
                return direction
            else:
                raise ValueError("Error, cannot work out gene's direction")
        else:
            return direction


class SiteFactory:
    """ A class used to represent a list of available sites (Site) on a DNA.

    Attributes
    ----------
    filename : str, optional
        Path to the site csv file.
    site_list : list
        A list containing Sites.
    """

    def __init__(self, filename=None):
        """ A class used to represent a list of available sites (Site) on a DNA.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file.
        """

        self.filename = filename
        self.site_list = []
        if filename:
            self.read_csv()

    def get_site_list(self):
        """ Gets the site_list

        Returns
        ----------
        list : A list of Sites.
        """
        return self.site_list

    def read_csv(self):
        """ Reads the SiteFactory csv filename and adds the sites to site_list.
        """
        df = pd.read_csv(self.filename)
        for index, row in df.iterrows():
            new_site = Site(site_type=row['type'], name=row['name'], start=float(row['start']), end=float(row['end']),
                            k_on=float(row['k_on']), binding_model_name=str(row['binding_model']),
                            binding_oparams_file=str(row['binding_oparams']))
            # new_site = Site(site_type=row['type'], name=row['name'], start=float(row['start']), end=float(row['end']),
            #  k_min=float(row['k_min']), k_max=float(row['k_max']), binding_model_name=row['model'],
            #                            binding_model_oparams=row['oparams'])
            self.site_list.append(new_site)
