import pandas as pd
from TORCphysics import binding_model as bm
from TORCphysics import effect_model as em
import sys
from TORCphysics import params


# TODO: Keep adding Effect models, then fix site, then enzyme. Remember, environment has binding, unbinding and enzyme
#  models. Enzymes have effect and unbinding. Sites only have binding models.
# TODO: If enviromentals are given, and no models are specified, the code should recommend/use default models
#  according the enzyme type

class Environment:
    """
    A class used to represent molecules/enzymes in the environment.

    Attributes
    ----------
    e_type : str
        The enzyme/molecule type, e.g. topo.
    name : str
        The name of the environmental/enzyme, e.g. gyrase
    site_list :
        list of sites that the enzyme can recognise and bind.
    concentration : float
        The concentration of the enzyme in nM.
    size : float
        The size of the enzyme in bp
    effective_size : float
        The effective size in bp. This size is assumed to be the size in which the enzyme makes contact with the DNA.
        So ef_size < size.
    site_type : str
    The type of site that this environmental can recognise to bind.
    binding_model_name : str, optional
         The name of the binding model to use
     binding_oparams_file : str, optional
         The path to the parameters for the binding model.
     effect_model_name : str, optional
         The name of the effect model to use
     effect_oparams_file : str, optional
         The path to the parameters for the effect model.
     unbinding_model_name : str, optional
         The name of the unbinding m odel to use
     unbinding_oparams_file : str, optional
         The path to the parameters for the ubinding model.
     binding_model : optional
         The binding model to use. This is a subclass of the BindingModel
     binding_model_oparams : dict, optional
         The actual dictionary with the parameters for the binding model
     effect_model : optional
         The effect model to use. This is a subclass of the EffectModel
     effect_model_oparams : dict, optional
         The actual dictionary with the parameters for the Effect model
     unbinding_model :optional
         The unbinding model to use. This is a subclass of the UnBindingModel
     unbinding_model_oparams : dict, optional
         The actual dictionary with the parameters for the unbinding model

    """

    def __init__(self, e_type, name, site_list, concentration, size, effective_size, site_type,
                 binding_model_name=None, binding_oparams_file=None,
                 effect_model_name=None, effect_oparams_file=None,
                 unbinding_model_name=None, unbinding_oparams_file=None,
                 binding_model=None, effect_model=None, unbinding_model=None,
                 binding_model_oparams=None, effect_model_oparams=None, unbinding_model_oparams=None):
        """
         A class used to represent molecules/enzymes in the environment.

         Attributes
         ----------
         e_type : str
             The enzyme/molecule type, e.g. topo.
         name : str
             The name of the environmental/enzyme, e.g. gyrase
         site_list :
             list of sites that the enzyme can recognise and bind.
         concentration : float
             The concentration of the enzyme in nM.
         size : float
             The size of the enzyme in bp
         effective_size : float
             The effective size in bp. This size is assumed to be the size in which the enzyme makes contact with the DNA.
             So ef_size < size.
         site_type : str
             The type of site that this environmental can recognise to bind.
         binding_model_name : str, optional
             The name of the binding model to use
         binding_oparams_file : str, optional
             The path to the parameters for the binding model.
         effect_model_name : str, optional
             The name of the effect model to use
         effect_oparams_file : str, optional
             The path to the parameters for the effect model.
         unbinding_model_name : str, optional
             The name of the unbinding m odel to use
         unbinding_oparams_file : str, optional
             The path to the parameters for the ubinding model.
         binding_model : optional
             The binding model to use. This is a subclass of the BindingModel
         binding_model_oparams : dict, optional
             The actual dictionary with the parameters for the binding model
         effect_model : optional
             The effect model to use. This is a subclass of the EffectModel
         effect_model_oparams : dict, optional
             The actual dictionary with the parameters for the Effect model
         unbinding_model :optional
             The unbinding model to use. This is a subclass of the UnBindingModel
         unbinding_model_oparams : dict, optional
             The actual dictionary with the parameters for the unbinding model
         """

        # Assign parameters
        print(name, binding_model_name, binding_oparams_file)
        self.enzyme_type = e_type
        self.name = name
        self.site_list = site_list  # It recognizes a list of sites, rather than a specific site
        self.site_type = site_type  # We need to remember the type
        self.concentration = concentration
        self.size = size
        self.effective_size = effective_size

        # Models
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

        self.check_inputs()

        self.get_models()

    #    def get_models(self, binding_model, effect_model, unbinding_model):

    def check_inputs(self):

        if not isinstance(self.enzyme_type, str) or self.enzyme_type == '':
            raise ValueError('Error, environmentals must have a type')
        if not isinstance(self.name, str) or self.name == '':
            raise ValueError('Error, environmentals must have a name')
        if not isinstance(self.site_list, list):
            raise ValueError('Error, environmentals site_list must be a list')
        if not isinstance(self.site_type, str) or self.site_type == '':
            raise ValueError('Error, environmentals need to recognise a site_type')
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
        if self.binding_model_name == '' or self.binding_model_name == 'None' or self.binding_model_name == 'none':
            self.binding_model_name = None
        if self.binding_model == '' or self.binding_model == 'None' or self.binding_model == 'none':
            self.binding_model = None
        if self.binding_oparams_file == '' or self.binding_oparams_file == 'None' or self.binding_oparams_file == 'none':
            self.binding_oparams_file = None
        if (self.binding_model_oparams == '' or self.binding_model_oparams == 'None' or
                self.binding_model_oparams == 'none'):
            self.binding_model_oparams = None

        # Effect model
        if self.effect_model_name == '' or self.effect_model_name == 'None' or self.effect_model_name == 'none':
            self.effect_model_name = None
        if self.effect_model == '' or self.effect_model == 'None' or self.effect_model == 'none':
            self.effect_model = None
        if self.effect_oparams_file == '' or self.effect_oparams_file == 'None' or self.effect_oparams_file == 'none':
            self.effect_oparams_file = None
        if (self.effect_model_oparams == '' or self.effect_model_oparams == 'None' or
                self.effect_model_oparams == 'none'):
            self.effect_model_oparams = None

        # Unbinding model
        if self.unbinding_model_name == '' or self.unbinding_model_name == 'None' or self.unbinding_model_name == 'none':
            self.unbinding_model_name = None
        if self.unbinding_model == '' or self.unbinding_model == 'None' or self.unbinding_model == 'none':
            self.unbinding_model = None
        if self.unbinding_oparams_file == '' or self.unbinding_oparams_file == 'None' or self.unbinding_oparams_file == 'none':
            self.unbinding_oparams_file = None
        if (self.unbinding_model_oparams == '' or self.unbinding_model_oparams == 'None' or
                self.unbinding_model_oparams == 'none'):
            self.unbinding_model_oparams = None

    # This function sorts the models
    def get_models(self):

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

            #  Let's check if it's actually a effect model - The model should already have the oparams
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
                    self.unbinding_model = bm.assign_unbinding_model(self.unbinding_model_name,
                                                                     **self.unbinding_model_oparams)
                # No dictionary was given
                else:
                    # If no oparams_file is given, then DEFAULT values are used.
                    if self.unbinding_oparams_file is None:
                        self.unbinding_model = bm.assign_unbinding_model(self.unbinding_model_name)
                    # If an oparams_file is given, then those are loaded
                    else:
                        self.unbinding_model = bm.assign_unbinding_model(self.unbinding_model_name,
                                                                         oparams_file=self.unbinding_oparams_file)

                    self.unbinding_model_oparams = self.unbinding_model.oparams  # To make them match

        # An actual model was given
        else:

            #  Let's check if it's actually a unbinding model - The model should already have the oparams
            if isinstance(self.unbinding_model, bm.UnBindingModel):
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
    def __init__(self, site_list, filename=None):
        self.filename = filename
        self.environment_list = []
        self.site_list = site_list
        if filename:
            self.read_csv()

    def get_environment_list(self):
        return self.environment_list

    def read_csv(self):
        df = pd.read_csv(self.filename)
        for index, row in df.iterrows():
            new_environment = Environment(e_type=row['type'], name=row['name'],
                                          site_list=self.site_match(row['site_type']),
                                          concentration=float(row['concentration']),
                                          size=float(row['size']),
                                          effective_size=float(row['effective_size']),
                                          site_type=row['site_type'],
                                          binding_model_name=row['binding_model'],
                                          binding_oparams_file=row['binding_oparams'],
                                          effect_model_name=row['effect_model'],
                                          effect_oparams_file=row['effect_oparams'],
                                          unbinding_model_name=row['unbinding_model'],
                                          unbinding_oparams_file=row['unbinding_oparams'])
            self.environment_list.append(new_environment)

    def site_match(self, label):
        #        enzyme_before = [enzyme.position for enzyme in enzyme_list if enzyme.position <= site.start][-1]
        site_list = [site for site in self.site_list if site.site_type == label]
        return site_list

#        if label in [site.name for site in self.site_list]:
#            for site in self.site_list:
#                if site.name == label:
#                    return site  # the first one?
#        else:
#            return None
