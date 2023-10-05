import numpy as np
from TORCphysics import params
import pandas as pd
from abc import ABC, abstractmethod


# ---------------------------------------------------------------------------------------------------------------------
# UNBINDING MODELS
# ---------------------------------------------------------------------------------------------------------------------
class UnBindingModel(ABC):
    """
     The UnBindingModel abstract class used for defining unbinding models (subclasses).
     If you need a new model, define it below.
     See how some of the models are defined from this class, so you can make your own and implement it.
    """

    def __init__(self, filename=None, **oparams):
        """ The constructor of UnBindingModel class.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the unbinding model.
        oparams : dict, optional
            A dictionary containing the parameters used for the unbinding model.
        """
        self.filename = filename
        self.oparams = oparams

    @abstractmethod
    def unbinding_probability(self) -> float:
        """ Abstract method for calculating the probability of unbinding.
        This is an essential function for UnBindingModels as they must be able to calculate the probability of
        unbinding for a given enzyme.

        Returns
        ----------
        probability : float
            It should return a probability (number), which indicates the probability of unbinding at the given timestep.
            Other functions/protocols should then interpret and implement this number.
        """
        pass


class PoissonUnBinding(UnBindingModel):
    """
     An unbinding model subclass that calculates unbinding probabilities according a Poisson process.
     This is one of the simplest unbinding models, where bound enzymes unbind at a constant rate.
    """

    def __init__(self, filename=None, **oparams):
        """ The constructor of the PoissonUnBinding subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the unbinding model; this file should have a k_off parameter
        oparams : dict, optional
            A dictionary containing the parameters used for the unbinding model. In this case, it would be k_off.
        """
        super().__init__(filename, **oparams)  # Call the base class constructor
        if not oparams:
            if filename is None:
                self.k_off = params.k_off
            else:
                rows = pd.read_csv(filename)
                self.k_off = float(rows['k_off'])
        else:
            self.k_off = oparams['k_off']

        self.oparams = {'k_off': self.k_off}  # Just in case

    #    def unbinding_probability(self, off_rate, dt) -> float:
    def unbinding_probability(self, dt) -> float:
        """ Method for calculating the probability of unbinding according a Poisson Process.

        Parameters
        ----------
        dt : float
            Timestep in seconds (s).

        Returns
        ----------
        probability : float
            A number that indicates the probability of unbinding in the current timestep.
        """
        return Poisson_process(self.k_off, dt)


# ---------------------------------------------------------------------------------------------------------------------
# HELPFUL FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------
# According inputs, loads the unbinding model, name and its params. This function is used in environment and enzymes.
# This function calls assign_unbinding_model
def get_unbinding_model(name, ub_model, model_name, oparams_file, oparams):
    """ This function loads the UnBindingModel to implement according the provided inputs.
    This function is used for environment and enzyme. So this function is implemented by those two classes.

    Parameters
    ----------
    name : str
        Name of the environmental or enzyme.
    ub_model : UnBindingModel or None
        A UnBindingModel or None.
    model_name : str
        Name of the model to use, e.g. 'PoissonUnBinding'
    oparams_file : str, optional
        Path to the csv file containing the parametrisation of the UnBindingModel to use.
    oparams : dict, optional
        A dictionary containing the parameters used for the unbinding model.
        In the case of PoissonUnBinding, it needs to have k_off.

    Returns
    ----------
    unbinding_model : UnBindingModel or None
        The UnBindingModel to implement for the Enzyme/Environment. If no UnBindingModel could be determined,
        this variable will be None.
    unbinding_model_name: str or None
        Name of the UnBindingModel to use. It is the same as unbinding_model.__class__.__name__
        If the UnBindingModel was not determined, then this variable is None.
    unbinding_oparams_file: str or None
        Path to the csv file containing the parametrisation of the UnBindingModel. None if file was not given.
    unbinding_model_oparams : dict or None
        Dictionary with the parametrisation of the UnBindingModel. None will be returned if the UnBindingModel could not
        be determined.
    """
    # If no model is given
    if ub_model is None:

        # No model is given, not even a name, so there's NO unbinding model
        if model_name is None:
            ub_model = None
            model_name = None
            oparams_file = None
            oparams = None

        # Model indicated by name
        else:
            # Loads unbinding model.
            # If oparams dict is given, those will be assigned to the model -> This is priority over oparams_file
            # If oparams_file is given, parameters will be read from file, in case of no oparams dict
            # If no oparams file/dict are given, default values will be used.

            # A dictionary of parameters is given so that's priority
            if isinstance(oparams, dict):
                ub_model = assign_unbinding_model(model_name, **oparams)
            # No dictionary was given
            else:
                # If no oparams_file is given, then DEFAULT values are used.
                if oparams_file is None:
                    ub_model = assign_unbinding_model(model_name)
                # If an oparams_file is given, then those are loaded
                else:
                    ub_model = assign_unbinding_model(model_name, oparams_file=oparams_file)

                oparams = ub_model.oparams  # To make them match

    # An actual model was given
    else:

        #  Let's check if it's actually an unbinding model - The model should already have the oparams
        if isinstance(ub_model, UnBindingModel):
            #  Then, some variables are fixed.
            model_name = ub_model.__class__.__name__
            oparams = ub_model.oparams
            oparams_file = None

        else:
            print('Warning, unbinding model given is not a class for environmental ', name)
            ub_model = None
            model_name = None
            oparams_file = None
            oparams = None

    unbinding_model = ub_model
    unbinding_model_name = model_name
    unbinding_oparams_file = oparams_file
    unbinding_model_oparams = oparams

    return unbinding_model, unbinding_model_name, unbinding_oparams_file, unbinding_model_oparams


# Add your models into this function so it the code can recognise it
def assign_unbinding_model(model_name, oparams_file=None, **oparams):
    """ This function decides the UnBindingModel to use according the provided inputs.

    Parameters
    ----------
    model_name : str
        Name of the UnBindingModel to use. e,g, PoissonUnBinding.
    oparams_file : str, optional
        Path to the csv file containing the parametrisation of the UnBindingModel to use.
    oparams : dict, optional
        A dictionary containing the parameters used for the unbinding model.
        In the case of PoissonUnBinding, it would need to have k_off.

    Returns
    ----------
    my_model : UnBindingModel
        A UnBindingModel object that describes the unbinding mechanism of the given site.
    """
    if model_name == 'PoissonUnBinding':
        my_model = PoissonUnBinding(filename=oparams_file, **oparams)
    else:
        raise ValueError('Could not recognise unbinding model ' + model_name)
    return my_model


def Poisson_process(rate, dt):
    rdt = rate * dt  # it is what is in the exponent (is that how you call it?)
    probability = rdt * np.exp(-rdt)
    return probability


def read_csv_to_dict(filename):
    """
    Reads csv file and puts it in a dictionary
    """
    return pd.read_csv(filename).to_dict()
