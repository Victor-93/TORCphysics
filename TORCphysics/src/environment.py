import pandas as pd
from TORCphysics import utils
from TORCphysics import binding_model as bm
from TORCphysics import unbinding_model as ubm
from TORCphysics import effect_model as em


class Environment:
    """
    A class used to represent molecules/enzymes in the environment.

    Attributes
    ----------
    enzyme_type : str
        The enzyme/molecule type, e.g. topo.
    name : str
        The name of the environmental/enzyme, e.g. gyrase
    site_list : list
        A list of sites that the enzyme can recognise and bind.
    concentration : float
        The concentration of the enzyme in nano Molars (nM).
    size : float
        The size of the enzyme in base-pairs (bp)
    effective_size : float
        The effective size in base-pairs (bp).
        This size is assumed to be the size in which the enzyme makes contact with the DNA.
        So, effective_size < size.
    site_type : str
        The type of site that this environmental can recognise and bind.
    binding_model_name : str, optional
        The name of the binding model to use, e.g., 'PoissonBinding'.
    binding_oparams_file : str, optional
        The path to the csv file containing the parameters for the binding model.
    binding_model : class, optional
        This is a subclass of the BindingModel to implement, e.g., PoissonBinding.
    binding_model_oparams : dict, optional
        A dictionary with the parameters for the binding model.
    effect_model_name : str, optional
        The name of the effect model to use, e.g., 'RNAPUniform'.
    effect_oparams_file : str, optional
        The path to the csv file containing the parameters for the effect model.
    effect_model : class, optional
        This is a subclass of EffectModel to implement, e.g., RNAPUniform.
    effect_model_oparams : dict, optional
        A dictionary with the parameters for the Effect model.
    unbinding_model_name : str, optional
        The name of the unbinding model to use, e.g., 'PoissonUnBinding'
    unbinding_oparams_file : str, optional
        The path to the csv file containing the parameters for the unbinding model.
    unbinding_model : class, optional
        This is a subclass of the UnBindingModel to use, e.g., PoissonUnBinding
    unbinding_model_oparams : dict, optional
        A actual dictionary with the parameters for the unbinding model.

    Notes
    ----------------
    Additional parameters (oparams) must be compatible with the given binding, effect or unbinding model, so please,
    do not create parameters that the models won't consider.
    When these oparams are not given, the code will load default parameters according the model indicated.
    These default parameters are saved in params.py
    """

    def __init__(self, e_type, name, site_list, concentration, size, effective_size, site_type,
                 binding_model_name=None, binding_oparams_file=None,
                 effect_model_name=None, effect_oparams_file=None,
                 unbinding_model_name=None, unbinding_oparams_file=None,
                 binding_model=None, effect_model=None, unbinding_model=None,
                 binding_model_oparams=None, effect_model_oparams=None, unbinding_model_oparams=None):
        """
        Constructor of the class Environment, used to represent the enzymes/molecules in the environment.

        Parameters
        ----------
        e_type : str
            The enzyme/molecule type, e.g. topo.
        name : str
            The name of the environmental/enzyme, e.g. gyrase
        site_list : list
            A list of sites that the enzyme can recognise and bind.
        concentration : float
            The concentration of the enzyme in nano Molars (nM).
        size : float
            The size of the enzyme in base-pairs (bp)
        effective_size : float
            The effective size in base-pairs (bp).
            This size is assumed to be the size in which the enzyme makes contact with the DNA.
            So, effective_size < size.
        site_type : str
            The type of site that this environmental can recognise and bind.
        binding_model_name : str, optional
            The name of the binding model to use, e.g., 'PoissonBinding'.
        binding_oparams_file : str, optional
            The path to the csv file containing the parameters for the binding model.
        binding_model : class, optional
            This is a subclass of the BindingModel to implement, e.g., PoissonBinding.
        binding_model_oparams : dict, optional
            A dictionary with the parameters for the binding model.
        effect_model_name : str, optional
            The name of the effect model to use, e.g., 'RNAPUniform'.
        effect_oparams_file : str, optional
            The path to the csv file containing the parameters for the effect model.
        effect_model : class, optional
            This is a subclass of EffectModel to implement, e.g., RNAPUniform.
        effect_model_oparams : dict, optional
            A dictionary with the parameters for the Effect model.
        unbinding_model_name : str, optional
            The name of the unbinding model to use, e.g., 'PoissonUnBinding'
        unbinding_oparams_file : str, optional
            The path to the csv file containing the parameters for the unbinding model.
        unbinding_model : class, optional
            This is a subclass of the UnBindingModel to use, e.g., PoissonUnBinding
        unbinding_model_oparams : dict, optional
            A actual dictionary with the parameters for the unbinding model.
         """

        # Assign parameters
        self.enzyme_type = e_type
        self.name = name
        self.site_list = site_list  # It recognizes a list of sites, rather than a specific site
        self.site_type = site_type  # We need to remember the type
        self.concentration = concentration
        self.size = size
        self.effective_size = effective_size

        # Assign models
        self.binding_model_name = binding_model_name
        self.binding_oparams_file = binding_oparams_file
        self.binding_model = binding_model
        self.binding_model_oparams = binding_model_oparams

        self.effect_model_name = effect_model_name
        self.effect_oparams_file = effect_oparams_file
        self.effect_model = effect_model
        self.effect_model_oparams = effect_model_oparams

        self.unbinding_model_name = unbinding_model_name
        self.unbinding_oparams_file = unbinding_oparams_file
        self.unbinding_model = unbinding_model
        self.unbinding_model_oparams = unbinding_model_oparams

        # Verify inputs
        self.check_inputs()

        # Loads the binding, effect and unbinding models if given.
        self.get_models()

    #    def get_models(self, binding_model, effect_model, unbinding_model):

    def check_inputs(self):
        """ Checks that Environment parameters are of the correct type.
        """

        if not isinstance(self.enzyme_type, str) or self.enzyme_type == '':
            raise ValueError('Error, environmentals must have a type')
        if not isinstance(self.name, str) or self.name == '':
            raise ValueError('Error, environmentals must have a name')
        if not isinstance(self.site_list, list):
            raise ValueError('Error, environmentals site_list must be a list')
        if (self.site_type == '' or self.site_type == 'None' or self.site_type == 'none' or
                self.site_type == 'nan'):
            self.site_type = ''
#        if not isinstance(self.site_type, str) or self.site_type == '':
#            raise ValueError('Error, environmentals need to recognise a site_type')
        if not isinstance(self.concentration, float) and not isinstance(self.concentration, int):
            raise ValueError('Error, environmentals need a number for concentration')
        if not isinstance(self.size, float) and not isinstance(self.size, int):
            raise ValueError('Error, environmentals need a number for size')
        if not isinstance(self.effective_size, float) and not isinstance(self.effective_size, int):
            raise ValueError('Error, environmentals need a number for effective size')
        if self.effective_size > self.size:
            print('Error, effective size effective_size cannot be larger than size')
            print('effective_size:', self.effective_size)
            print('size', self.size)
            print('For environmental ', self.name)
            raise ValueError('Error: effective_size > size')

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

        # Effect model
        if (self.effect_model_name == '' or self.effect_model_name == 'None' or self.effect_model_name == 'none' or
                self.effect_model_name == 'nan'):
            self.effect_model_name = None
        if (self.effect_model == '' or self.effect_model == 'None' or self.effect_model == 'none' or
                self.effect_model == 'nan'):
            self.effect_model = None
        if (self.effect_oparams_file == '' or self.effect_oparams_file == 'None' or self.effect_oparams_file == 'none'
                or self.effect_oparams_file == 'nan'):
            self.effect_oparams_file = None
        if (self.effect_model_oparams == '' or self.effect_model_oparams == 'None' or
                self.effect_model_oparams == 'none' or self.effect_model_oparams == 'nan'):
            self.effect_model_oparams = None

        # Unbinding model
        if (self.unbinding_model_name == '' or self.unbinding_model_name == 'None'
                or self.unbinding_model_name == 'none' or self.unbinding_model_name == 'nan'):
            self.unbinding_model_name = None
        if (self.unbinding_model == '' or self.unbinding_model == 'None' or self.unbinding_model == 'none' or
                self.unbinding_model == 'nan'):
            self.unbinding_model = None
        if (self.unbinding_oparams_file == '' or self.unbinding_oparams_file == 'None'
                or self.unbinding_oparams_file == 'none' or self.unbinding_oparams_file == 'nan'):
            self.unbinding_oparams_file = None
        if (self.unbinding_model_oparams == '' or self.unbinding_model_oparams == 'None' or
                self.unbinding_model_oparams == 'none' or self.unbinding_model_oparams == 'nan'):
            self.unbinding_model_oparams = None

    def get_models(self):
        """ Loads the Environment's binding, effect and unbinding models (if given).
        """

        # Binding Model
        self.binding_model, self.binding_model_name, self.binding_oparams_file, self.binding_model_oparams = (
            bm.get_binding_model(self.name, self.binding_model, self.binding_model_name,
                                 self.binding_oparams_file, self.binding_model_oparams))
        # Effect Model
        self.effect_model, self.effect_model_name, self.effect_oparams_file, self.effect_model_oparams = (
            em.get_effect_model(self.name, self.effect_model, self.effect_model_name,
                                self.effect_oparams_file, self.effect_model_oparams))

        # Unbinding Model
        self.unbinding_model, self.unbinding_model_name, self.unbinding_oparams_file, self.unbinding_model_oparams = (
            ubm.get_unbinding_model(self.name, self.unbinding_model, self.unbinding_model_name,
                                    self.unbinding_oparams_file, self.unbinding_model_oparams))

    # PREVIOUS VERSION. IT IS OBSOLETE NOW, SO PLEASE REMOVE IT
    # This function sorts the models
    def get_models2(self):

        # Binding model
        # -------------------------------------------------------------
        # If no model is given
        if self.binding_model is None:

            # No model is given, not even a name, so there's NO binding model
            if self.binding_model_name is None:
                self.binding_model = None
                self.binding_model_name = None
                self.binding_oparams_file = None
                self.binding_model_oparams = None

            # Model indicated by name
            else:
                # Loads binding model.
                # If oparams dict is given, those will be assigned to the model -> This is priority over oparams_file
                # If oparams_file is given, parameters will be read from file, in case of no oparams dict
                # If no oparams file/dict are given, default values will be used.

                # A dictionary of parameters is given so that's priority
                if isinstance(self.binding_model_oparams, dict):
                    self.binding_model = bm.assign_binding_model(self.binding_model_name,
                                                                 **self.binding_model_oparams)
                # No dictionary was given
                else:
                    # If no oparams_file is given, then DEFAULT values are used.
                    if self.binding_oparams_file is None:
                        self.binding_model = bm.assign_binding_model(self.binding_model_name)
                    # If an oparams_file is given, then those are loaded
                    else:
                        self.binding_model = bm.assign_binding_model(self.binding_model_name,
                                                                     oparams_file=self.binding_oparams_file)

                    self.binding_model_oparams = self.binding_model.oparams  # To make them match

        # An actual model was given
        else:

            #  Let's check if it's actually a binding model - The model should already have the oparams
            if isinstance(self.binding_model, bm.BindingModel):
                #  Then, some variables are fixed.
                self.binding_model_name = self.binding_model.__class__.__name__
                self.binding_model_oparams = self.binding_model.oparams
                self.binding_oparams_file = None

            else:
                print('Warning, binding model given is not a class for environmental ', self.name)
                self.binding_model = None
                self.binding_model_name = None
                self.binding_oparams_file = None
                self.binding_model_oparams = None

        # Effect model
        # -------------------------------------------------------------
        # If no model is given
        if self.effect_model is None:

            # No model is given, not even a name, so there's NO effect model
            if self.effect_model_name is None:
                self.effect_model = None
                self.effect_model_name = None
                self.effect_oparams_file = None
                self.effect_model_oparams = None

            # Model indicated by name
            else:
                # Loads effect model.
                # If oparams dict is given, those will be assigned to the model -> This is priority over oparams_file
                # If oparams_file is given, parameters will be read from file, in case of no oparams dict
                # If no oparams file/dict are given, default values will be used.

                # A dictionary of parameters is given so that's priority
                if isinstance(self.effect_model_oparams, dict):
                    self.effect_model = em.assign_effect_model(self.effect_model_name,
                                                               **self.effect_model_oparams)
                # No dictionary was given
                else:
                    # If no oparams_file is given, then DEFAULT values are used.
                    if self.effect_oparams_file is None:
                        self.effect_model = em.assign_effect_model(self.effect_model_name)
                    # If an oparams_file is given, then those are loaded
                    else:
                        self.effect_model = em.assign_effect_model(self.effect_model_name,
                                                                   oparams_file=self.effect_oparams_file)

                    self.effect_model_oparams = self.effect_model.oparams  # To make them match

        # An actual model was given
        else:

            #  Let's check if it's actually an effect model - The model should already have the oparams
            if isinstance(self.effect_model, em.EffectModel):
                #  Then, some variables are fixed.
                self.effect_model_name = self.effect_model.__class__.__name__
                self.effect_model_oparams = self.effect_model.oparams
                self.effect_oparams_file = None

            else:
                print('Warning, effect model given is not a class for environmental ', self.name)
                self.effect_model = None
                self.effect_model_name = None
                self.effect_oparams_file = None
                self.effect_model_oparams = None

        # Unbinding model
        # -------------------------------------------------------------
        # If no model is given
        if self.unbinding_model is None:

            # No model is given, not even a name, so there's NO unbinding model
            if self.unbinding_model_name is None:
                self.unbinding_model = None
                self.unbinding_model_name = None
                self.unbinding_oparams_file = None
                self.unbinding_model_oparams = None

            # Model indicated by name
            else:
                # Loads unbinding model.
                # If oparams dict is given, those will be assigned to the model -> This is priority over oparams_file
                # If oparams_file is given, parameters will be read from file, in case of no oparams dict
                # If no oparams file/dict are given, default values will be used.

                # A dictionary of parameters is given so that's priority
                if isinstance(self.unbinding_model_oparams, dict):
                    self.unbinding_model = ubm.assign_unbinding_model(self.unbinding_model_name,
                                                                      **self.unbinding_model_oparams)
                # No dictionary was given
                else:
                    # If no oparams_file is given, then DEFAULT values are used.
                    if self.unbinding_oparams_file is None:
                        self.unbinding_model = ubm.assign_unbinding_model(self.unbinding_model_name)
                    # If an oparams_file is given, then those are loaded
                    else:
                        self.unbinding_model = ubm.assign_unbinding_model(self.unbinding_model_name,
                                                                          oparams_file=self.unbinding_oparams_file)

                    self.unbinding_model_oparams = self.unbinding_model.oparams  # To make them match

        # An actual model was given
        else:

            #  Let's check if it's actually an unbinding model - The model should already have the oparams
            if isinstance(self.unbinding_model, ubm.UnBindingModel):
                #  Then, some variables are fixed.
                self.unbinding_model_name = self.unbinding_model.__class__.__name__
                self.unbinding_model_oparams = self.unbinding_model.oparams
                self.unbinding_oparams_file = None

            else:
                print('Warning, unbinding model given is not a class for environmental ', self.name)
                self.unbinding_model = None
                self.unbinding_model_name = None
                self.unbinding_oparams_file = None
                self.unbinding_model_oparams = None


class EnvironmentFactory:
    """ A class used to represent a list of environmentals (enzymes/molecules) in the environment.

    Attributes
    ----------
    site_list : list
        A list containing Sites on the DNA.
    environment_list : list
        A list containing Environmentals, or in other words, the molecules/enzymes on the environment.
    filename : str, optional
        Path to the environment csv file.
    """

    def __init__(self, site_list, filename=None):
        """ Constructor of the class EnvironmentFactory.

        Parameters
        ----------
        site_list : list
            A list containing Sites on the DNA.
        filename : str, optional
            Path to the environment csv file.
        """
        self.filename = filename
        self.environment_list = []
        self.site_list = site_list
        if filename:
            self.read_csv()

    def get_environment_list(self):
        """ Gets the environment_list

        Returns
        ----------
        list : A list of environmentals (enzymes/molecules) in the environment.
        """

        return self.environment_list

    def read_csv(self):
        """ Reads the EnvironmentFactory csv filename and adds the environmentals to environment_list.
        """
        df = pd.read_csv(self.filename)
        for index, row in df.iterrows():
            new_environment = Environment(e_type=str(row['type']), name=str(row['name']),
                                          site_list=utils.site_match_by_type(site_list=self.site_list,
                                                                             label=row['site_type']),
                                          concentration=float(row['concentration']),
                                          size=float(row['size']),
                                          effective_size=float(row['effective_size']),
                                          site_type=str(row['site_type']),
                                          binding_model_name=str(row['binding_model']),
                                          binding_oparams_file=str(row['binding_oparams']),
                                          effect_model_name=str(row['effect_model']),
                                          effect_oparams_file=str(row['effect_oparams']),
                                          unbinding_model_name=str(row['unbinding_model']),
                                          unbinding_oparams_file=str(row['unbinding_oparams']))
            self.environment_list.append(new_environment)
