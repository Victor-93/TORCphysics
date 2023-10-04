import sys
import numpy as np
from TORCphysics import params
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
# UNBINDING MODELS
# ---------------------------------------------------------------------------------------------------------------------
# These functions are for applying particular unbinding conditions.
# Their output is "unbind", which can be True if the bound enzyme will unbind, or False if it will not unbind (stays).
# If you need to add your own, create it here!

# Model that controls the unbinding condition for DNA topoisomerases
# ---------------------------------------------------------------------------------------------------------------------
def RNAP_unbinding_model(enzyme):
    # condition for transcription in >>>>>>>>>>>>> right direction or
    # condition for transcription in <<<<<<<<<<<<< left  direction
    unbind = False
    if (enzyme.direction == 1 and enzyme.end - enzyme.position <= 0) or \
            (enzyme.direction == -1 and enzyme.end - enzyme.position >= 0):
        unbind = True
    return unbind


# Unbinding model for enzymes that unbind spontaneously
# ---------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------
# BINDING MODELS
# ---------------------------------------------------------------------------------------------------------------------
# According inputs, loads the binding model, name and its params. This function is used in environment and sites.
# This function calls assign_binding_model
def get_binding_model(name, b_model, model_name, oparams_file, oparams):
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
            print('Warning, binding model given is not a class for environmental ', name)
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
    if model_name == 'PoissonBinding':
        my_model = PoissonBinding(filename=oparams_file, **oparams)
    elif model_name == 'TopoIRecognition':
        my_model = TopoIRecognition(filename=oparams_file, **oparams)
    #    elif model_name == 'GyraseRecognition':
    #        my_model = GyraseRecognition(filename=oparams_file)
    else:
        raise ValueError('Could not recognise binding model ' + model_name)
    return my_model


# These functions are for the particular binding model.
# If you need a new one, define it here.
class BindingModel(ABC):
    """
     The BindingModel abstract class used for defining binding models (subclasses).
     If you need a new model, define it below.
     See how some of the models are defined from this class, so you can make your own and implement it.
    """

    def __init__(self, filename=None, **oparams):
        """ The constructor of BindingModel.

        Parameters
        ----------
        filename : str, optional
            Path to the site csv file that parametrises the binding model.
        oparams : dict, optional
            A dictionary containing the parameters used for the binding model.
        """

        self.filename = filename
        self.oparams = oparams

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


# Model for binding probability according Poisson process
class PoissonBinding(BindingModel):
    """
     A binding model subclass that calculates binding probabilities according a Poisson process.
     This is one of the simplest binding models, where enzymes bind at a constant rate.
    """

    def __init__(self, filename=None, **oparams):
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
                rows = pd.read_csv(filename)
                self.k_on = float(rows['k_on'])
        else:
            self.k_on = oparams['k_on']

        self.oparams = {'k_on': self.k_on}  # Just in case

    def binding_probability(self, on_rate, dt) -> float:
        """ Method for calculating the probability of binding according a Poisson Process.

        Returns
        ----------
        probability : float
            A number that indicates the probability of binding in the current timestep.
        """

        return Poisson_process(on_rate, dt)


class TopoIRecognition(BindingModel):
    def __init__(self, filename=None, **oparams):
        super().__init__(filename, **oparams)  # name  # Call the base class constructor
        if not oparams:
            if filename is None:
                self.width = params.topo_b_w
                self.threshold = params.topo_b_t
                self.k_on = params.topo_b_k_on
            else:
                rows = pd.read_csv(filename)
                self.width = float(rows['width'])
                self.threshold = float(rows['threshold'])
                self.k_on = float(rows['k_on'])
        else:
            self.width = float(oparams['width'])
            self.threshold = float(oparams['threshold'])
            self.k_on = float(oparams['k_on'])

        self.oparams = {'width': self.width, 'threshold': self.threshold, 'k_on': self.k_on}  # Just in case

    # Notice that the concentration of enzyme is outside the model as it can vary during the simulation.
    def binding_probability(self, enzyme, sigma) -> float:

        a = enzyme.concentration * self.k_on
        b = 1 + np.exp((sigma - self.threshold) / self.width)
        rate = a / b
        return rate


class GyraseRecognition(BindingModel):

    def __init__(self, filename=None, **oparams):
        super().__init__(filename, **oparams)  # name  # Call the base class constructor
        if not oparams:
            if filename is None:
                self.width = params.gyra_b_w
                self.threshold = params.gyra_b_t
                self.k_on = params.gyra_b_k_on
            else:
                rows = pd.read_csv(filename)
                self.width = float(rows['width'])
                self.threshold = float(rows['threshold'])
                self.k_on = float(rows['k_on'])
        else:
            self.width = float(oparams['width'])
            self.threshold = float(oparams['threshold'])
            self.k_on = float(oparams['k_on'])

        self.oparams = {'width': self.width, 'threshold': self.threshold, 'k_on': self.k_on}  # Just in case

    # Notice that the concentration of enzyme is outside the model as it can vary during the simulation.
    def binding_probability(self, enzyme, sigma) -> float:

        a = enzyme.concentration * self.k_on
        b = 1 + np.exp(-(sigma - self.threshold) / self.width)
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


# ----------------------------------------------------------
# This equation calculates the probability of binding according
# a Non-homogeneous Poisson process, which is basically a Poisson process
# with variable rate (simply modelling).
# It assumes that the rate was already varied and obtained by one of the opening energies
# sigma - supercoiling density
def P_binding_Nonh_Poisson(rate, dt):
    probability = rate * dt  # The smaller dt the more accurate it is.

    return probability


# ----------------------------------------------------------
# The promoter activation curve according Sam Meyer 2019
# For this function, we use the minimum rate
def promoter_curve_Meyer(basal_rate, sigma):
    u = 1.0 / (1.0 + np.exp((sigma - SM_sigma_t) / SM_epsilon_t))  # the energy required for melting
    f = np.exp(SM_m * u)  # the activation curve
    rate = basal_rate * f  # and the rate modulated through the activation curve
    return rate


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


def Poisson_process(rate, dt):
    rdt = rate * dt  # it is what is in the exponent (is that how you call it?)
    probability = rdt * np.exp(-rdt)
    return probability


def read_csv_to_dict(filename):
    """
    Reads csv file and puts it in a dictionary
    """
    return pd.read_csv(filename).to_dict()


# This function checks if a site is not blocked by other enzymes
def check_site_availability(site, enzyme_list, size):
    # Check if the site is available for binding.
    # It assumes that the probability has already been calculated, and we have a candidate enzyme for the binding
    # with size=size.
    # We need the list of current enzymes to see if the one before and after the site overlap with the start site.
    enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= site.start][-1]
    enzyme_after = [enzyme for enzyme in enzyme_list if enzyme.position >= site.start][0]
    # And a range of their occupancy
    range_before = [enzyme_before.position, enzyme_before.position + enzyme_before.size]
    range_after = [enzyme_after.position, enzyme_after.position + enzyme_after.size]
    # TODO: Check if this is correct! I think it is assuming that enzymes bind just before the start site, which might
    #  not be true.
    if site.direction > 0:
        my_range = [site.start - size, site.start]
    #        my_range = [site.start, site.start + size]
    elif site.direction <= 0:
        my_range = [site.start, site.start + size]
    else:
        print('Error in checking site availability. Site=', site.site_type, site.name)
        sys.exit()
        #  my_range = [site.start, site.start + size]
    #        my_range = [site.start, site.start - size]

    # If any of them intersect
    if (set(range_before) & set(my_range)) or (set(range_after) & set(my_range)):
        available = False
    # there is an intersection
    else:
        available = True

    return available


# ----------------------------------------------------------
# This function makes sure that only one enzyme will end up binding a region.
# It checks that the enzymes in the list of new_enzymes do not overlap and if they do, decide which will end up
# binding
# TODO: test this function - design a experiment in which you kind of know what outcome you should get.
def check_binding_conflicts(enzyme_list, rng):
    enzyme_list.sort(key=lambda x: x.position)  # sort by position
    checked_enzyme_list = []
    s = 0
    for i, my_enzyme in enumerate(enzyme_list):
        if i == 0:  # We need enzymes after
            continue
        enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= my_enzyme.position][-1]

        # Check if they overlap
        if enzyme_before.position + enzyme_before.size >= my_enzyme.position:
            # It overlaps, so decide which enzymes stays
            urandom = rng.uniform()  # we need a random number for the decision
            if urandom <= 0.5:  # If it is <= 0.5, then the enzyme before stays.
                del enzyme_list[i - s - 1]
                s += 1
                # checked_enzyme_list.append(enzyme_before)
            else:  # And if >0.5, then we don't add the enzyme before (we lose it).
                del enzyme_list[i - s]
                s += 1
                # continue
        # else:
        # checked_enzyme_list.append(enzyme_before)  # If nothing overlaps, then nothing happens

    return enzyme_list  # checked_enzyme_list
