import pandas as pd
from TORCphysics import binding_model as bm
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
    eff_size : float
        The effective size in bp. This size is assumed to be the size in which the enzyme makes contact with the DNA.
        So ef_size < size.
    site_type : str
        The type of site that this environmental can recognize to bind.
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
    effect_model : optional
        The effect model to use. This is a subclass of the EffectModel
    unbinding_model :optional
        The unbinding model to use. This is a subclass of the UnBindingModel

    Methods
    -------
    says(sound=None)
        Prints the animals name and what sound it makes
    """
    def __init__(self, e_type, name, site_list, concentration, size, eff_size, site_type,
                 binding_model_name=None, binding_oparams_file=None,
                 effect_model_name=None, effect_oparams_file=None,
                 unbinding_model_name=None, unbinding_oparams_file=None,
                 binding_model=None, effect_model=None, unbinding_model=None):
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
         eff_size : float
             The effective size in bp. This size is assumed to be the size in which the enzyme makes contact with the DNA.
             So ef_size < size.
         site_type : str
             The type of site that this environmental can recognize to bind.
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
         effect_model : optional
             The effect model to use. This is a subclass of the EffectModel
         unbinding_model :optional
             The unbinding model to use. This is a subclass of the UnBindingModel
         """

        #        self.binding_oparams = None
        #        self.effect_oparams = None
        # TODO: Finish this class, then test it! YOU NEED TO COMPLETELY CHECK IT
        # TODO: Add the option where you select default models? Or noup? Maybe no for the moment.
        self.binding_model = None
        self.effect_model = None
        self.unbinding_model = None
        self.enzyme_type = e_type
        self.name = name
        self.site_list = site_list  # It recognizes a list of sites, rather than a specific site
        self.site_type = site_type  # We need to remember the type
        self.concentration = concentration
        self.size = size
        self.binding_model_name = binding_model_name
        self.binding_oparams_file = binding_oparams_file
        self.effect_model_name = effect_model_name
        self.effect_oparams_file = effect_oparams_file
        self.unbinding_model_name = unbinding_model_name
        self.unbinding_oparams_file = unbinding_oparams_file

    def get_models(self, binding_model, effect_model, unbinding_model):

        # Binding model
        if issubclass(binding_model, bm.BindingModel):
            self.binding_model = binding_model
        if ((binding_model is None)
                and (self.binding_model_name != 'none'
                     or self.binding_model_name != 'None'
                     or self.binding_model_name is not None)):
            self.binding_model = bm.assign_binding_model(self.binding_model_name)

        # Effect Model

        # Unbinding Model
        if issubclass(unbinding_model, bm.UnBindingModel):
            self.unbinding_model = unbinding_model
        if ((unbinding_model is None)
                and (self.unbinding_model_name != 'none'
                     or self.unbinding_model_name != 'None'
                     or self.unbinding_model_name is not None)):
            self.unbinding_model = bm.assign_unbinding_model(self.unbinding_model_name)


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
                                          concentration=float(row['concentration']), size=float(row['size']),
                                          site_type=row['site_type'],
                                          binding_model=row['binding_model'],
                                          binding_oparams=row['binding_oparams'],
                                          effect_model=row['effect_model'],
                                          effect_oparams=row['effect_oparams'],
                                          unbinding_model=row['unbinding_model'],
                                          unbinding_oparams=row['unbinding_oparams'])
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
