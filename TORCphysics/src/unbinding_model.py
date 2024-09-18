import numpy as np
from TORCphysics import params, utils
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

     Attributes
     ----------
     filename : str, optional
         Path to the site csv file that parametrises the unbinding model.
     oparams : dict, optional
         A dictionary containing the parameters used for the unbinding model.
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
     An UnbindingModel subclass that calculates unbinding probabilities according a Poisson process.
     This is one of the simplest unbinding models, where bound enzymes unbind at a constant rate.

     Attributes
     ----------
     k_off : float
        Rate (1/s) at which the enzymes unbind.
     filename : str, optional
         Path to the site csv file that parametrises the unbinding model.
     oparams : dict, optional
         A dictionary containing the parameters used for the unbinding model.
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
                mydata = pd.read_csv(filename)
                if 'k_off' in mydata.columns:
                    #  self.k_on = float(rows['k_on'])
                    self.k_off = mydata['k_off'][0]
                else:
                    raise ValueError('Error, k_off parameter missing in csv file for PoissonUnBinding')
        else:
            self.k_off = oparams['k_off']

        self.oparams = {'k_off': self.k_off}  # Just in case

    #    def unbinding_probability(self, off_rate, dt) -> float:
    def unbinding_probability(self, enzyme, dt) -> float:
        """ Method for calculating the probability of unbinding according a Poisson Process.

        Parameters
        ----------
        enzyme : Enzyme
            The enzyme that is currently transcribing the DNA and its site.
        dt : float
            Timestep in seconds (s).

        Returns
        ----------
        probability : float
            A number that indicates the probability of unbinding in the current timestep.
        """
        return utils.Poisson_process(self.k_off, dt)

# TODO: Gyrase unbinding that depends and takes into account the energy.


class RNAPSimpleUnbinding(UnBindingModel):
    """
     An UnbindingModel subclass that calculates unbinding probabilities of a RNAP transcribing a gene.
     In this simple model, the RNAP unbinds if it reaches the end of the transcribing region.

     This model doesn't really have any useful attributes.

     Attributes
     ----------
     filename : str, optional
         Path to the site csv file that parametrises the unbinding model.
     oparams : dict, optional
         A dictionary containing the parameters used for the unbinding model.
    """

    def __init__(self, filename=None, **oparams):
        """ The constructor of the PoissonUnBinding subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the unbinding model; for this model, there's nothing useful in
            a file.
        oparams : dict, optional
            A dictionary containing the parameters used for the unbinding model. In this case, it is not required
        """
        super().__init__(filename, **oparams)  # Call the base class constructor

    #    def unbinding_probability(self, off_rate, dt) -> float:
    def unbinding_probability(self, enzyme, dt) -> float:
        """ Method for calculating the probability of unbinding according a Poisson Process.

        Parameters
        ----------
        enzyme : Enzyme
            The enzyme that is currently transcribing the DNA and its site.
        dt : float
            Timestep in seconds (s).

        Returns
        ----------
        probability : float
            A number that indicates the probability of unbinding in the current timestep.
        """

        probability = 0.0  # It will not unbind unless fulfils the condition

        if enzyme.direction != 1 and enzyme.direction != -1:
            raise ValueError('Error. Enzyme with invalid direction in RNAPSimpleUnbinding.')

        # TODO: Check how it changes with the sizes and effective_size
        # condition for transcription in >>>>>>>>>>>>> right direction or
        # condition for transcription in <<<<<<<<<<<<< left  direction
        if (enzyme.direction == 1 and enzyme.end - enzyme.position <= 0) or \
                (enzyme.direction == -1 and enzyme.end - enzyme.position >= 0):
            probability = 1.0  # There's no way it won't unbind with this

        return probability

class RNAPStagesSimpleUnbinding(UnBindingModel):
    """
     An UnbindingModel subclass that calculates unbinding probabilities of a RNAP bound to the DNA.
     There are two possible scenarios in this model in which unbinding can happen.
     1. If the enzyme is in the Closed_complex state, then it can spontaniusly unbind the DNA according a Poisson
     process.
     2. If the enzyme is transcribing the DNA and reaches the end of the transcribing region, it will immediately
     unbind the DNA.

     Attributes
     ----------
     k_off : float
        Rate (1/s) at which the enzymes unbind when in the Closed_complex state.
     filename : str, optional
         Path to the site csv file that parametrises the unbinding model.
     oparams : dict, optional
         A dictionary containing the parameters used for the unbinding model.
    """

    def __init__(self, filename=None, **oparams):
        """ The constructor of the PoissonUnBinding subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the unbinding model; for this model, there's nothing useful in
            a file.
        oparams : dict, optional
            A dictionary containing the parameters used for the unbinding model. In this case, it is not required
        """
        super().__init__(filename, **oparams)  # Call the base class constructor
        if not oparams:
            if filename is None:
                self.k_off = params.k_off
            else:
                mydata = pd.read_csv(filename)
                if 'k_off' in mydata.columns:
                    #  self.k_on = float(rows['k_on'])
                    self.k_off = mydata['k_off'][0]
                else:
                    raise ValueError('Error, k_off parameter missing in csv file for RNAPStagesSimpleUnbinding.')
        else:
            self.k_off = oparams['k_off']

        self.oparams = {'k_off': self.k_off}  # Just in case


    #    def unbinding_probability(self, off_rate, dt) -> float:
    def unbinding_probability(self, enzyme, dt) -> float:
        """ Method for calculating the probability of unbinding according a Poisson Process.

        Parameters
        ----------
        enzyme : Enzyme
            The enzyme that is currently transcribing the DNA and its site.
        dt : float
            Timestep in seconds (s).

        Returns
        ----------
        probability : float
            A number that indicates the probability of unbinding in the current timestep.
        """

        probability = 0.0  # It will not unbind unless fulfils the condition

        # If in the closed complex state, it can spontaneously unbind
        if enzyme.effect_model.state == 'Closed_complex':
            probability = utils.Poisson_process(self.k_off, dt)

        # If it is in the elongation state (transcribing), it'll only unbind if it reaches the terminator/end.
        if enzyme.effect_model.state == 'Elongation':

            if enzyme.direction != 1 and enzyme.direction != -1:
                raise ValueError('Error. Enzyme with invalid direction in RNAPStagesSimpleUnbinding.')

            # condition for transcription in >>>>>>>>>>>>> right direction or
            # condition for transcription in <<<<<<<<<<<<< left  direction
            if (enzyme.direction == 1 and enzyme.end - enzyme.position <= 0) or \
                    (enzyme.direction == -1 and enzyme.end - enzyme.position >= 0):
                probability = 1.0  # There's no way it won't unbind with this

        return probability

# This one should replace the one above. It is essentially the same model/function, but it communicates with the
# site to get the k_off.
class RNAPStagesSimpleUnbindingv2(UnBindingModel):
    """
     An UnbindingModel subclass that calculates unbinding probabilities of a RNAP bound to the DNA.
     There are two possible scenarios in this model in which unbinding can happen.
     1. If the enzyme is in the Closed_complex state, then it can spontaniusly unbind the DNA according a Poisson
     process.
     2. If the enzyme is transcribing the DNA and reaches the end of the transcribing region, it will immediately
     unbind the DNA.

     Attributes
     ----------
     filename : str, optional
         Path to the site csv file that parametrises the unbinding model.
     oparams : dict, optional
         A dictionary containing the parameters used for the unbinding model.
    """

    def __init__(self, filename=None, **oparams):
        """ The constructor of the PoissonUnBinding subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the unbinding model; for this model, there's nothing useful in
            a file.
        oparams : dict, optional
            A dictionary containing the parameters used for the unbinding model. In this case, it is not required
        """
        super().__init__(filename, **oparams)  # Call the base class constructor

        # It doesn't have any related parameters on its own


    #    def unbinding_probability(self, off_rate, dt) -> float:
    def unbinding_probability(self, enzyme, dt) -> float:
        """ Method for calculating the probability of unbinding according a Poisson Process.

        Parameters
        ----------
        enzyme : Enzyme
            The enzyme that is currently transcribing the DNA and its site.
        dt : float
            Timestep in seconds (s).

        Returns
        ----------
        probability : float
            A number that indicates the probability of unbinding in the current timestep.
        """

        probability = 0.0  # It will not unbind unless fulfils the condition
        k_off = enzyme.site.binding_model.k_off

        # If in the closed complex state, it can spontaneously unbind
        if enzyme.effect_model.state == 'Closed_complex':
            probability = utils.Poisson_process(k_off, dt)

        # If it is in the elongation state (transcribing), it'll only unbind if it reaches the terminator/end.
        if enzyme.effect_model.state == 'Elongation':

            if enzyme.direction != 1 and enzyme.direction != -1:
                raise ValueError('Error. Enzyme with invalid direction in RNAPStagesSimpleUnbinding.')

            # condition for transcription in >>>>>>>>>>>>> right direction or
            # condition for transcription in <<<<<<<<<<<<< left  direction
            if (enzyme.direction == 1 and enzyme.end - enzyme.position <= 0) or \
                    (enzyme.direction == -1 and enzyme.end - enzyme.position >= 0):
                probability = 1.0  # There's no way it won't unbind with this

        return probability

# TODO: Document and test
class LacISimpleUnBinding(UnBindingModel):

    def __init__(self, filename=None, **oparams):

        super().__init__(filename, **oparams)  # Call the base class constructor
        if not oparams:
            if filename is None:
                self.k_off = 0.01  # Rate at which enzymes unbind
            else:
                mydata = pd.read_csv(filename)
                if 'k_off' in mydata.columns:
                    self.k_off = mydata['k_off'][0]
                else:
                    raise ValueError('Error, k_off parameter missing in csv file for LacISimpleUnBinding.')
        else:
            self.k_off = oparams['k_off']

        self.oparams = {'k_off': self.k_off}  # Just in case

    def unbinding_probability(self, enzyme, dt) -> float:

        # If the bridge is formed, then it can't unbind
        if enzyme.effect_model.state == 'UNLOOPED' :   #NOT JUST_UNLOOPED, because we remain with one
            probability = utils.Poisson_process(self.k_off, dt)
        # If the bridge isn't formed, then it can unbind
        else:
            probability = 0.0
        return probability



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
    elif model_name == 'RNAPSimpleUnbinding':
        my_model = RNAPSimpleUnbinding(filename=oparams_file, **oparams)
    elif model_name == 'RNAPStagesSimpleUnbinding':
        my_model = RNAPStagesSimpleUnbinding(filename=oparams_file, **oparams)
    elif model_name == 'RNAPStagesSimpleUnbindingv2':
        my_model = RNAPStagesSimpleUnbindingv2(filename=oparams_file, **oparams)
    elif model_name == 'LacISimpleUnBinding':
        my_model = LacISimpleUnBinding(filename=oparams_file, **oparams)
    else:
        raise ValueError('Could not recognise unbinding model ' + model_name)
    return my_model


