import numpy as np
from TORCphysics import params
import pandas as pd
from abc import ABC, abstractmethod


# ---------------------------------------------------------------------------------------------------------------------
# UNBINDING MODELS
# ---------------------------------------------------------------------------------------------------------------------
# According inputs, loads the unbinding model, name and its params. This function is used in environment and sites.
# This function calls assign_unbinding_model
def get_unbinding_model(name, ub_model, model_name, oparams_file, oparams):
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
    if model_name == 'PoissonUnBinding':
        my_model = PoissonUnBinding(filename=oparams_file, **oparams)
    else:
        raise ValueError('Could not recognise unbinding model ' + model_name)
    return my_model


class UnBindingModel(ABC):
    def __init__(self, filename=None, **oparams):
        self.filename = filename
        self.oparams = oparams

    @abstractmethod
    def unbinding_probability(self) -> float:
        pass


class PoissonUnBinding(UnBindingModel):
    def __init__(self, filename=None, **oparams):
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

    def unbinding_probability(self, off_rate, dt) -> float:
        return Poisson_process(off_rate, dt)


# ---------------------------------------------------------------------------------------------------------------------
# HELPFUL FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------

def Poisson_process(rate, dt):
    rdt = rate * dt  # it is what is in the exponent (is that how you call it?)
    probability = rdt * np.exp(-rdt)
    return probability


def read_csv_to_dict(filename):
    """
    Reads csv file and puts it in a dictionary
    """
    return pd.read_csv(filename).to_dict()