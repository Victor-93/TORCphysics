import sys
import numpy as np
from TORCphysics import params, utils
import pandas as pd
from abc import ABC, abstractmethod

# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This module contains the mathematical functions that compose the statistical part of my model
# It comprises the necessary equations required for simulating the stochastic binding, and the
# topoisomerases/gyrase activities

# All parameters are already in the params module, but I prefer to have them here with more simple names:
v0 = params.v0
w0 = params.w0
gamma = params.gamma
# dt     = params.dt

kBT = 310.0 * params.kB_kcalmolK  # The Boltzmann constant multiplied by 310K which is the temperature
# at which the SIDD code is ran...

# Sam Meyer's PROMOTER CURVE (parameters taken from Houdaigi NAR 2019)
SM_sigma_t = params.SM_sigma_t
SM_epsilon_t = params.SM_epsilon_t
SM_m = params.SM_m

# EFFECTIVE ENERGY PROMOTER CURVE (inspired by Houdaigi NAR 2019)
EE_alpha = params.EE_alpha

# Topoisomerase activity parameters
topo_w = params.topo_b_w
topo_t = params.topo_b_t
gyra_w = params.gyra_b_w
gyra_t = params.gyra_b_t


# TODO: We still need to check if the new enzymes are not overlapping. If more than 1 enzyme
#  passed the probability test and are binding the same overlapping region, we need to flip a coin to
#  check which one will be binding


# ---------------------------------------------------------------------------------------------------------------------
# BINDING MODELS
# ---------------------------------------------------------------------------------------------------------------------
# TODO: I added the parameter interacts. Think about it and see if you remove it or maybe not?
# These functions are for the particular binding model.
# If you need a new one, define it here.
class BindingModel(ABC):
    """
     The BindingModel abstract class used for defining binding models (subclasses).
     If you need a new model, define it below.
     See how some of the models are defined from this class, so you can make your own and implement it.

     Attributes
     ----------
     filename : str, optional
         Path to the site csv file that parametrises the binding model.
     interacts : bool, optional
        Parameter that indicates if the site/environmental interacts with other bound enzymes. It could be the case
        that the bound enzymem may increase the rate of the binding environmental.
     oparams : dict, optional
         A dictionary containing the parameters used for the binding model.
    """

    def __init__(self, filename=None, interacts=False, **oparams):
        """ The constructor of BindingModel.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the binding model.
        interacts : bool, optional
            Parameter that indicates if the site/environmental interacts with other bound enzymes. It could be the case
            that the bound enzymem may increase the rate of the binding environmental.
        oparams : dict, optional
            A dictionary containing the parameters used for the binding model.
        """

        self.filename = filename
        self.oparams = oparams
        self.interacts = interacts

    @abstractmethod
    def binding_probability(self) -> float:
        """ Abstract method for calculating the probability of binding.
        This is an essential function for BindingModels as they must be able to calculate the probability of binding of
        a given enzyme.

        Returns
        ----------
        probability : float
            It should return a probability (number), which indicates the probability of binding at the given timestep.
            Other functions/protocols should then interpret and implement this number.
        """
        pass

    # TODO: Document this abstract method
    @abstractmethod
    def rate_modulation(self) -> float:
        pass


# TODO: Maybe we need to add the reference concentration! OR remember that the concentration you give, you always
#       need to divide it by the reference in which the model was calibrated.
#       Or maybe not if you calibrate doing k_on * concentration
# Model for binding probability according Poisson process
class PoissonBinding(BindingModel):
    """
     A BindingModel subclass that calculates binding probabilities according a Poisson process.
     This is one of the simplest binding models, where enzymes bind at a constant rate.

     Attributes
     ----------
     k_on : float
        Rate (1/s) at which the enzymes bind.
     filename : str, optional
        Path to the site csv file that parametrises the binding model.
     oparams : dict, optional
        A dictionary containing the parameters used for the binding model.
    """

    def __init__(self, filename=None, interacts=False, **oparams):
        """ The constructor of the PoissonBinding subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the binding model; this file should have a k_on parameter
        oparams : dict, optional
            A dictionary containing the parameters used for the binding model. In this case, it would be k_on.
        """
        super().__init__(filename, **oparams)  # Call the base class constructor
        if not oparams:
            if filename is None:
                self.k_on = params.k_on  # If no filename was given, then it uses default k_on.
                # Note that the value used also depends on the input site.csv file.
            else:
                mydata = pd.read_csv(filename)
                if 'k_on' in mydata.columns:
                    #  self.k_on = float(rows['k_on'])
                    self.k_on = mydata['k_on'][0]
                else:
                    raise ValueError('Error, k_on parameter missing in csv file for PoissonBinding')
        else:
            self.k_on = oparams['k_on']

        self.interacts = interacts
        self.oparams = {'k_on': self.k_on}  # Just in case

    #    def binding_probability(self, on_rate, dt) -> float:
    # NOTE: Shouldn't on_rate be the same that k_on? It should be a property of the Site on it's own, right?
    def binding_probability(self, environmental, superhelical, dt) -> float:
        """ Method for calculating the probability of binding according a Poisson Process.

        Parameters
        ----------
        dt : float
            Timestep in seconds (s).
        environmental : Environment
            The environmental molecule that is trying to bind the site.
        superhelical : float
            The local supercoiling region in the site.

        Returns
        ----------
        probability : float
            A number that indicates the probability of binding in the current timestep.
        """
        return utils.Poisson_process(self.k_on * environmental.concentration, dt)

    def rate_modulation(self, superhelical) -> float:
        rate = self.k_on
        return rate


# TODO: Document this function
class MeyerPromoterOpening(BindingModel):

    def __init__(self, filename=None, interacts=False, **oparams):

        super().__init__(filename, **oparams)  # Call the base class constructor
        if not oparams:
            if filename is None:
                self.k_on = params.k_on  # If no filename was given, then it uses default k_on.
                # Note that the value used also depends on the input site.csv file.
            else:
                mydata = pd.read_csv(filename)
                if 'k_on' in mydata.columns:
                    #  self.k_on = float(rows['k_on'])
                    self.k_on = mydata['k_on'][0]
                else:
                    raise ValueError('Error, k_on parameter missing in csv file for MeyerPromoterOpening')
        else:
            self.k_on = oparams['k_on']

        self.interacts = interacts
        self.oparams = {'k_on': self.k_on}  # Just in case

    #    def binding_probability(self, on_rate, dt) -> float:
    # NOTE: Shouldn't on_rate be the same that k_on? It should be a property of the Site on it's own, right?
    def binding_probability(self, environmental, superhelical, dt) -> float:
        rate = utils.promoter_curve_Meyer(basal_rate=self.k_on, superhelical=superhelical)
        probability = utils.P_binding_Nonh_Poisson(rate=rate, dt=dt)
        return probability

    def rate_modulation(self, superhelical) -> float:
        rate = utils.promoter_curve_Meyer(basal_rate=self.k_on, superhelical=superhelical)
        return rate


# TODO: Do you need the inverse version of the function? And if so, would you need another class or is it better
#  with a boolean?
#  Also check overflows/underflows
# TODO: Needs testing
class MaxMinPromoterBinding(BindingModel):
    """
     A BindingModel subclass that calculates binding probabilities according a sigmoid response curve.
     This model represents the binding mechanism of supercoiling sensitive promoters, which openning curve usually
     follows a sigmoidal response as function of the local supercoiling density. In this case, the binding rate varies
     as a function of supercoiling

     Attributes
     ----------
     k_min : float
        Minimum rate (1/s) at which the enzymes bind.
     k_max : float
        Maximum rate (1/s) at which the enzymes bind.
     threshold : float
        The threshold of the sigmoid curve. This is a dimensionless parameter.
     width : float
        The width of the sigmoid curve. This is a dimensionless parameter.
     filename : str, optional
        Path to the site csv file that parametrises the binding model.
     oparams : dict, optional
        A dictionary containing the parameters used for the binding model.
    """

    def __init__(self, filename=None, interacts=False, **oparams):

        super().__init__(filename, **oparams)
        if not oparams:
            if filename is None:
                self.width = .003
                self.threshold = -.094
                self.k_min = 0.001
                self.k_max = 0.01
            else:
                mydata = pd.read_csv(filename)
                if 'k_min' in mydata.columns:
                    self.k_min = mydata['k_min'][0]
                else:
                    raise ValueError('Error, k_max parameter missing in csv file for MaxMinPromoterBinding')
                if 'k_max' in mydata.columns:
                    self.k_max = mydata['k_max'][0]
                else:
                    raise ValueError('Error, k_max parameter missing in csv file for MaxMinPromoterBinding')
                if 'width' in mydata.columns:
                    self.width = mydata['width'][0]
                else:
                    raise ValueError('Error, width parameter missing in csv file for MaxMinPromoterBinding')
                if 'threshold' in mydata.columns:
                    self.threshold = mydata['threshold'][0]
                else:
                    raise ValueError('Error, threshold parameter missing in csv file for MaxMinPromoterBinding')
        else:
            # No point testing or checking that we have these variables, as python gives error on its own.
            self.width = float(oparams['width'])
            self.threshold = float(oparams['threshold'])
            self.k_min = float(oparams['k_min'])
            self.k_max = float(oparams['k_max'])

        self.interacts = interacts
        # Just in case
        self.oparams = {'width': self.width, 'threshold': self.threshold, 'k_min': self.k_min, 'k_max': self.k_max}

        # Notice that the concentration of enzyme is outside the model as it can vary during the simulation.

    def binding_probability(self, environmental, superhelical, dt) -> float:
        """ Method for calculating the probability of binding according the MaxMinPromoter model.

        Parameters
        ----------
        dt : float
            Timestep in seconds (s).
        environmental : Environment
            The environmental molecule that is trying to bind the site.
        superhelical : float
            The local supercoiling region in the site.

        Returns
        ----------
        probability : float
            A number that indicates the probability of binding in the current timestep.
        """

        a = np.log(self.k_min / self.k_max)
        b = 1 + np.exp(-(superhelical - self.threshold) / self.width)
        rate = self.k_max * np.exp(a / b)
        return utils.P_binding_Nonh_Poisson(rate=rate, dt=dt)

    def rate_modulation(self, superhelical) -> float:
        """ Method for calculating the rate modulation as a functino of superhelical density for the
        MaxMinPromoter model.

        Parameters
        ----------
        dt : float
            Timestep in seconds (s).
        environmental : Environment
            The environmental molecule that is trying to bind the site.
        superhelical : float
            The local supercoiling region in the site.

        Returns
        ----------
        rate : float
            Rate at which the environmentals bind given the local superhelical density
        """

        a = np.log(self.k_min / self.k_max)
        b = 1 + np.exp(-(superhelical - self.threshold) / self.width)
        rate = self.k_max * np.exp(a / b)
        return rate


class TopoIRecognition(BindingModel):
    """
     A BindingModel subclass that calculates binding probabilities according a sigmoid recognition curve.
     This model represents the binding mechanism of Topoisomerase I, where it recognises the DNA's shape according
     the local supercoiling density. In this case, the binding rate varies as a function of supercoiling

     Attributes
     ----------
     k_on : float
        Rate (1/s) at which the enzymes bind.
     threshold : float
        The threshold of the sigmoid curve. This is a dimensionless parameter.
     width : float
        The width of the sigmoid curve. This is a dimensionless parameter.
     filename : str, optional
        Path to the site csv file that parametrises the binding model.
     oparams : dict, optional
        A dictionary containing the parameters used for the binding model.
    """

    def __init__(self, filename=None, interacts=False, **oparams):
        """ The constructor of the TopoIRecognition subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the binding model; this file should have the k_on, width and
             threshold parameters.
        oparams : dict, optional
            A dictionary containing the parameters used for the binding model. In this case, it would be k_on, width
            and threshold.
        """

        super().__init__(filename, **oparams)
        if not oparams:
            if filename is None:
                self.width = params.topo_b_w
                self.threshold = params.topo_b_t
                self.k_on = params.topo_b_k_on
            else:
                mydata = pd.read_csv(filename)
                if 'k_on' in mydata.columns:
                    self.k_on = mydata['k_on'][0]
                else:
                    raise ValueError('Error, k_on parameter missing in csv file for TopoIRecognition')
                if 'width' in mydata.columns:
                    self.width = mydata['width'][0]
                else:
                    raise ValueError('Error, width parameter missing in csv file for TopoIRecognition')
                if 'threshold' in mydata.columns:
                    self.threshold = mydata['threshold'][0]
                else:
                    raise ValueError('Error, threshold parameter missing in csv file for TopoIRecognition')
        else:
            # No point testing or checking that we have these variables, as python gives error on its own.
            self.width = float(oparams['width'])
            self.threshold = float(oparams['threshold'])
            self.k_on = float(oparams['k_on'])

        self.interacts = interacts
        self.oparams = {'width': self.width, 'threshold': self.threshold, 'k_on': self.k_on}  # Just in case

    # Notice that the concentration of enzyme is outside the model as it can vary during the simulation.
    def binding_probability(self, environmental, superhelical, dt) -> float:
        """ Method for calculating the probability of binding according the TopoIRecognition model.

        Parameters
        ----------
        dt : float
            Timestep in seconds (s).
        environmental : Environment
            The environmental molecule that is trying to bind the site.
        superhelical : float
            The local supercoiling region in the site.

        Returns
        ----------
        probability : float
            A number that indicates the probability of binding in the current timestep.
        """

        a = environmental.concentration * self.k_on
        #b = 1 + np.exp((superhelical - self.threshold) / self.width)
        #b = 1 + np.exp((superhelical + self.threshold) / self.width)
        b1 = np.float128((superhelical + self.threshold) / self.width)
        b = 1 + np.exp(b1)
        rate = a / b
        return utils.P_binding_Nonh_Poisson(rate=rate, dt=dt)

    def rate_modulation(self, superhelical) -> float:
        # Note that is not multiplied by the concentration
        a = self.k_on
        #b = 1 + np.exp((superhelical - self.threshold) / self.width)
        b = 1 + np.exp((superhelical + self.threshold) / self.width)
        rate = a / b
        return rate


class TopoIRecognitionRNAPTracking(BindingModel):
    """
     A BindingModel subclass that calculates binding probabilities according a sigmoid recognition curve and a
     RNAP tracking function.
     This model represents the binding mechanism of Topoisomerase I, where it recognises the DNA's shape according
     the local supercoiling density and the presence of RNAPs facilitates the binding (increases basal rate).
     In this case, the binding rate varies as a function of supercoiling and RNAP position.

     Attributes
     ----------
     k_on : float
        Rate (1/s) at which the enzymes bind.
     threshold : float
        The threshold of the sigmoid curve. This is a dimensionless parameter.
     width : float
        The width of the sigmoid curve. This is a dimensionless parameter.
    RNAP_dist: float
        Distance (bp) at which bound RNAPs affect the binding of topo I.
    fold_change: float
        Fold change due to the presence of RNAPs. Basically, if a topo I binds within a region in which a RNAP
        is located within the "RNAP_dist" distance, then the rate is increased by: rate x fold_change.
     filename : str, optional
        Path to the site csv file that parametrises the binding model.
     oparams : dict, optional
        A dictionary containing the parameters used for the binding model.
    """

    def __init__(self, filename=None, interacts=True, **oparams):
        """ The constructor of the TopoIRecognition subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the binding model; this file should have the k_on, width and
             threshold parameters.
        oparams : dict, optional
            A dictionary containing the parameters used for the binding model. In this case, it would be k_on, width
            and threshold.
        """

        super().__init__(filename, interacts, **oparams)
        if not oparams:
            if filename is None:
                self.width = params.topo_b_w
                self.threshold = params.topo_b_t
                self.k_on = params.topo_b_k_on
            else:
                mydata = pd.read_csv(filename)
                if 'k_on' in mydata.columns:
                    self.k_on = mydata['k_on'][0]
                else:
                    raise ValueError('Error, k_on parameter missing in csv file for TopoIRecognition')
                if 'width' in mydata.columns:
                    self.width = mydata['width'][0]
                else:
                    raise ValueError('Error, width parameter missing in csv file for TopoIRecognition')
                if 'threshold' in mydata.columns:
                    self.threshold = mydata['threshold'][0]
                else:
                    raise ValueError('Error, threshold parameter missing in csv file for TopoIRecognition')
        else:
            # No point testing or checking that we have these variables, as python gives error on its own.
            self.width = float(oparams['width'])
            self.threshold = float(oparams['threshold'])
            self.k_on = float(oparams['k_on'])

        self.RNAP_dist = 200.0  # in bp
        self.fold_change = 10#40.0  # Fold change when RNAP is close
        self.interacts = interacts
        self.oparams = {'width': self.width, 'threshold': self.threshold, 'k_on': self.k_on}  # Just in case

    # Notice that the concentration of enzyme is outside the model as it can vary during the simulation.
    def binding_probability(self, environmental, superhelical, site, enzyme_list, dt) -> float:
        """ Method for calculating the probability of binding according the TopoIRecognitionRNAPTracking model.

        Parameters
        ----------
        dt : float
            Timestep in seconds (s).
        environmental : Environment
            The environmental molecule that is trying to bind the site.
        superhelical : float
            The local supercoiling region in the site.

        Returns
        ----------
        probability : float
            A number that indicates the probability of binding in the current timestep.
        """

        # Calculate rate based on superhelical density
        a = environmental.concentration * self.k_on
        b1 = np.float128((superhelical + self.threshold) / self.width)
        b = 1 + np.exp(b1)
        rate = a / b
        # b = 1 + np.exp((superhelical - self.threshold) / self.width)


        # See if RNAP is in the vicinity and if yes, then the rate increases.
        # TODO:
        #   1.- We need to check that the RNAP is a distance less or equal than RNAP_dist.
        #   2.- Then check that has the correct direction:
        #    * If left and direction == -1, then increases
        #    * If right and direction == +1, then increases
        #  Necesitamos distancias del site:  s1____s2
        #  Necesitamos
        RNAP_near = False
        if site.start < site.end:
            s1 = site.start
            s2 = site.end
        else:
            s1 = site.end
            s2 = site.start

        # For Enzyme on the left
        enzyme_left = utils.get_enzyme_before_position(s1, enzyme_list)
        if 'RNAP' in enzyme_left.enzyme_type:
            dist_left = s1 - (enzyme_left.position + enzyme_left.effective_size)
            if 0 <= dist_left <= self.RNAP_dist and enzyme_left.direction == -1:
                RNAP_near = True

        # For Enzyme on the right
        enzyme_right = utils.get_enzyme_after_position(s2, enzyme_list)
        if 'RNAP' in enzyme_right.enzyme_type:
            dist_right = enzyme_right.position - s2
            if 0 <= dist_right <= self.RNAP_dist and enzyme_right.direction == 1:
                RNAP_near = True

        # Check if rate is increased
        if RNAP_near:
            rate = rate * self.fold_change
            #print('yes', probability *self.fold_change)
        # else:
            #print('no', probability)
        probability = utils.P_binding_Nonh_Poisson(rate=rate, dt=dt)
        return probability

    def rate_modulation(self, superhelical) -> float:
        # Note that is not multiplied by the concentration
        a = self.k_on
        b = 1 + np.exp((superhelical - self.threshold) / self.width)
        rate = a / b
        return rate


class GyraseRecognition(BindingModel):
    """
     A BindingModel subclass that calculates binding probabilities according a sigmoid recognition curve.
     This model represents the binding mechanism of Gyrase, where it recognises the DNA's shape according
     the local supercoiling density. In this case, the binding rate varies as a function of supercoiling

     Attributes
     ----------
     k_on : float
        Rate (1/s) at which the enzymes bind.
     threshold : float
        The threshold of the sigmoid curve. This is a dimensionless parameter.
     width : float
        The width of the sigmoid curve. This is a dimensionless parameter.
     filename : str, optional
        Path to the site csv file that parametrises the binding model.
     oparams : dict, optional
        A dictionary containing the parameters used for the binding model.
    """

    def __init__(self, filename=None, interacts=False, **oparams):
        """ The constructor of the GyraseRecognition subclass.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the binding model; this file should have the k_on, width and
             threshold parameters.
        oparams : dict, optional
            A dictionary containing the parameters used for the binding model. In this case, it would be k_on, width
            and threshold.
        """

        super().__init__(filename, **oparams)
        if not oparams:
            if filename is None:
                self.width = params.gyra_b_w
                self.threshold = params.gyra_b_t
                self.k_on = params.gyra_b_k_on
            else:
                mydata = pd.read_csv(filename)
                if 'k_on' in mydata.columns:
                    self.k_on = mydata['k_on'][0]
                else:
                    raise ValueError('Error, k_on parameter missing in csv file for GyraseRecognition')
                if 'width' in mydata.columns:
                    self.width = mydata['width'][0]
                else:
                    raise ValueError('Error, width parameter missing in csv file for GyraseRecognition')
                if 'threshold' in mydata.columns:
                    self.threshold = mydata['threshold'][0]
                else:
                    raise ValueError('Error, threshold parameter missing in csv file for GyraseRecognition')
        else:
            self.width = float(oparams['width'])
            self.threshold = float(oparams['threshold'])
            self.k_on = float(oparams['k_on'])

        self.interacts = interacts
        self.oparams = {'width': self.width, 'threshold': self.threshold, 'k_on': self.k_on}  # Just in case

    # Notice that the concentration of enzyme is outside the model as it can vary during the simulation.
    def binding_probability(self, environmental, superhelical, dt) -> float:
        """ Method for calculating the probability of binding according the GyraseRecognition model.

        Parameters
        ----------
        dt : float
            Timestep in seconds (s).
        environmental : Environment
            The environmental molecule that is trying to bind the site.
        superhelical : float
            The local supercoiling region in the site.

        Returns
        ----------
        probability : float
            A number that indicates the probability of binding in the current timestep.
        """

        a = environmental.concentration * self.k_on
        b1 = np.float128(-(superhelical + self.threshold) / self.width)
        b = 1 + np.exp(b1)
        # b = 1 + np.exp(-(superhelical - self.threshold) / self.width)
        #b = 1 + np.exp(-(superhelical + self.threshold) / self.width)
        #sigma_cons = .5
        #b = 1 + np.exp(-(superhelical + self.threshold) / self.width - sigma_cons) * np.exp(sigma_cons)
        rate = a / b
        # Configure NumPy to suppress overflow warnings
        #np.seterr(over='ignore')

        #try:
         #   # Your original expression
         #   a = environmental.concentration * self.k_on
         #   b1 = np.exp(-(superhelical + self.threshold) / self.width)
        #    # b = 1 +
            #rate = a / b
        #except:
        #    # Handle overflow condition
        #    print('Overflow in GyraseRecognition, superhelical = ', superhelical)
        #    rate = 0  # Or any other strategy you prefer
        #    sys.exit()
        #else:
        #    b = 1 + b1
        #    rate = a/b
        return utils.P_binding_Nonh_Poisson(rate=rate, dt=dt)

    def rate_modulation(self, superhelical) -> float:
        # Note that is not multiplied by the concentration
        a = self.k_on
        # b = 1 + np.exp(-(superhelical - self.threshold) / self.width)
        b = 1 + np.exp(-(superhelical + self.threshold) / self.width)
        rate = a / b
        return rate


# class lacI_binding(binding_model):
#    def __init__(self, name, bridging, unbridging):
#        super().__init__(name)  # Call the base class constructor
#        self.bridging = bridging
#        self.unbridging = unbridging

# Creating instances of the subclasses
# topoI_instance = topoI_binding(name="TopoI Binding", width=10, threshold=5)
# lacI_instance = lacI_binding(name="lacI Binding", bridging=True, unbridging=False)

# ---------------------------------------------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------
# According inputs, loads the binding model, name and its params. This function is used in environment and sites.
# This function calls assign_binding_model
def get_binding_model(name, b_model, model_name, oparams_file, oparams):
    """ This function loads the BindingModel to implement according the provided inputs.
    This function is used for environment and sites. So this function is implemented by those two classes.

    Parameters
    ----------
    name : str
        Name of the environmental or site.
    b_model : BindingModel or None
        A BindingModel or None.
    model_name : str
        Name of the model to use, e.g. 'PoissonBinding'
    oparams_file : str, optional
        Path to the csv file containing the parametrisation of the BindingModel to use.
    oparams : dict, optional
        A dictionary containing the parameters used for the binding model. In this case, it would be k_on.

    Returns
    ----------
    binding_model : BindingModel or None
        The BindingModel to implement for the Site/Environment. If no BindingModel could be determined, this variable
        will be None.
    binding_model_name: str or None
        Name of the BindingModel to use. It is the same as binding_model.__class__.__name__
        If the BindingModel was not determined, then this variable is None.
    binding_oparams_file: str or None
        Path to the csv file containing the parametrisation of the BindingModel. None if file was not given.
    binding_model_oparams : dict or None
        Dictionary with the parametrisation of the BindingModel. None will be returned if the BindingModel could not
        be determined.
    """

    # If no model is given
    if b_model is None:

        # No model is given, not even a name, so there's NO binding model
        if model_name is None:
            b_model = None
            model_name = None
            oparams_file = None
            oparams = None

        # Model indicated by name
        else:
            # Loads binding model.
            # If oparams dict is given, those will be assigned to the model -> This is priority over oparams_file
            # If oparams_file is given, parameters will be read from file, in case of no oparams dict
            # If no oparams file/dict are given, default values will be used.

            # A dictionary of parameters is given so that's priority
            if isinstance(oparams, dict):
                b_model = assign_binding_model(model_name, **oparams)
            # No dictionary was given
            else:
                # If no oparams_file is given, then DEFAULT values are used.
                if oparams_file is None:
                    b_model = assign_binding_model(model_name)
                # If an oparams_file is given, then those are loaded
                else:
                    b_model = assign_binding_model(model_name, oparams_file=oparams_file)

                oparams = b_model.oparams  # To make them match

    # An actual model was given
    else:

        #  Let's check if it's actually a binding model - The model should already have the oparams
        if isinstance(b_model, BindingModel):
            #  Then, some variables are fixed.
            model_name = b_model.__class__.__name__
            oparams = b_model.oparams
            oparams_file = None

        else:
            print('Warning, binding model given is not a class for environmental or site ', name)
            b_model = None
            model_name = None
            oparams_file = None
            oparams = None

    binding_model = b_model
    binding_model_name = model_name
    binding_oparams_file = oparams_file
    binding_model_oparams = oparams

    return binding_model, binding_model_name, binding_oparams_file, binding_model_oparams


# Add your models into this function so it the code can recognise it
def assign_binding_model(model_name, oparams_file=None, **oparams):
    """ This function decides the BindingModel to use according the provided inputs.

    Parameters
    ----------

    model_name : str
        Name of the BindingModel to use. e,g, PoissonBinding.
    oparams_file : str, optional
        Path to the csv file containing the parametrisation of the BindingModel to use.
    oparams : dict, optional
        A dictionary containing the parameters used for the binding model. In this case, it would be k_on.

    Returns
    ----------
    my_model : BindingModel
        A BindingModel object that describes the binding mechanism of the given site.
    """

    if model_name == 'PoissonBinding':
        my_model = PoissonBinding(filename=oparams_file, **oparams)
    elif model_name == 'TopoIRecognition':
        my_model = TopoIRecognition(filename=oparams_file, **oparams)
    elif model_name == 'TopoIRecognitionRNAPTracking':
        my_model = TopoIRecognitionRNAPTracking(filename=oparams_file, **oparams)
    elif model_name == 'GyraseRecognition':
        my_model = GyraseRecognition(filename=oparams_file, **oparams)
    elif model_name == 'MeyerPromoterOpening':
        my_model = MeyerPromoterOpening(filename=oparams_file, **oparams)
    elif model_name == 'MaxMinPromoterBinding':
        my_model = MaxMinPromoterBinding(filename=oparams_file, **oparams)
    else:
        raise ValueError('Could not recognise binding model ' + model_name)
    return my_model


# ----------------------------------------------------------


# ----------------------------------------------------------
# The supercoiling dependant opening energy of the promoter
# sequence. It follows a sigmoid curve, and should be
# calculated with the SIDD algorithm and following the methods
# of Houdaigui 2021 for isolating the discriminator sequence.

# Parameters:
# sigma - supercoiling density
# a,b,sigma_t,epsilon - sigmoid curve fitted parameters
def opening_energy(x, a, b, sigma_t, epsilon):
    return a + b / (1 + np.exp(-(x - sigma_t) / epsilon))


# ----------------------------------------------------------
# The promoter activation curve relaying on the effective
# thermal energy. This curve is parametrized by the fitting
# of the openning energy, and inspired by according Sam Meyer 2019
# For this function, we use the minimum rate
def promoter_curve_opening_E(basal_rate, sigma, sigma0, *opening_p):
    u = opening_energy(sigma, *opening_p)  # the energy required for melting
    u0 = opening_energy(sigma0, *opening_p)  # energy for melting at reference sigma0
    # (should be the sigma at which k0=basal_rate
    #  was measured...)
    du = u - u0  # Energy difference
    f = np.exp(-du / EE_alpha)  # the activation curve
    rate = basal_rate * f
    return rate


# ----------------------------------------------------------
# The promoter activation curve parametrized by the
# opening energy fitted parameters, and the observed
# maximum and minimum rates.
# opening_p - opening energy fitted parameters
# k_min = minimum rate
# k_max = maximum rate
def promoter_curve_opening_E_maxmin(k_min, k_max, sigma, *opening_p):
    a = np.log(k_min / k_max)
    b = 1 + np.exp(-(sigma - opening_p[2]) / opening_p[3])
    rate = k_max * np.exp(a / b)
    return rate


# ----------------------------------------------------------
# Basically, is the same function as the previous one,
# but this one has the sigmoid inverted, hence, we
# need to adjust the location of the minimum/maximum rates in
# the equation.
# opening_p - opening energy fitted parameters
# k_min = minimum rate
# k_max = maximum rate
def promoter_curve_opening_E_maxmin_I(k_min, k_max, sigma, *opening_p):
    a = np.log(k_max / k_min)
    b = 1 + np.exp(-(sigma - opening_p[2]) / opening_p[3])
    rate = k_min * np.exp(a / b)
    return rate


# TODO: Check the curve doesn't overflow?
def topoI_binding(topo, sigma):
    a = topo.concentration * topo.k_on
    # TODO: Check why sometimes the len is bigger
    b = 1 + np.exp((sigma - topo.oparams['threshold']) / topo.oparams['width'])  # [0] because the dictionary
    #    b = 1 + np.exp((sigma - topo.oparams['threshold'][0]) / topo.oparams['width'][0]) # [0] because the dictionary
    # gives me trouble
    rate = a / b  # * np.exp(1 / b)
    #    try:
    #        b = 1 + np.exp((sigma - topo_t) / topo_w)
    #        rate = a / b
    #    except OverflowError as oe:
    #        sigma_removed = 0.0
    return rate


def gyrase_binding(gyra, sigma):
    a = gyra.concentration * gyra.k_on
    b = 1 + np.exp(-(sigma - gyra.oparams['threshold']) / gyra.oparams['width'])
    #    b = 1 + np.exp(-(sigma - gyra.oparams['threshold'][0]) / gyra.oparams['width'][0]) # [0] because the dictionary
    # gives me trouble
    #    rate = a * np.exp(1 / b)
    rate = a / b
    return rate

# ---------------------------------------------------------------------------------------------------------------------
# HELPFUL FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------
# These functions are used as auxiliary.
# Can be for checking site availability, if multiple enzymes are trying to bind the same site, decide which one
# will bind, etc...

# ----------------------------------------------------------
