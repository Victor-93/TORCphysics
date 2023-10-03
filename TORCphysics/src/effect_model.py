import numpy as np
from TORCphysics import params
import pandas as pd
from abc import ABC, abstractmethod
import sys

# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This module contains the mathematical functions that compose the effects model.
# Effects can describe the motion of RNAPs as well as twist injection; topoisomerase activity; overall mechanics of
# enzymes bound to DNA, etc...

# TODO: Decide which of these parameters you need
# All paramete# rs are already in the params module, but I prefer to have them here with more simple names:
v0 = params.v0
w0 = params.w0
gamma = params.gamma


# ---------------------------------------------------------------------------------------------------------------------
# EFFECT
# ---------------------------------------------------------------------------------------------------------------------
class Effect:
    # I thought it'll be easier to describe the effects as an object.
    # Because the current effects are taken place at the current enzyme i=index, I use the index to locate the enzyme
    # in the enzyme_list which is applying the effect.
    # These effects act locally, so they can modify the enzyme's position, and the twist at the neighbouring domains.
    def __init__(self, index, position, twist_left, twist_right):
        # I'll save the input filenames just in case
        self.index = index
        self.position = position
        self.twist_left = twist_left
        self.twist_right = twist_right


# ---------------------------------------------------------------------------------------------------------------------
# EFFECT MODELS
# ---------------------------------------------------------------------------------------------------------------------
class EffectModel(ABC):
    """
     The EffectModel class for defining effect models.
     If you need a new one, define it here.
    """

    def __init__(self, filename=None, **oparams):
        self.filename = filename
        self.oparams = oparams

    @abstractmethod
    def calculate_effect(self) -> Effect:
        """
        Effects models need a "calculate_effect" function. This function must return an Effect object
        """
        pass


# ----------------------
# YOU can define your own models here!
# ----------------------

class RNAPUniform(EffectModel):
    """
     Effect model for the RNAP uniform motion. RNAP moves with constant speed and injects +-supercoils on each site
    """

    # def __init__(self, name, filename):
    def __init__(self, filename=None, **oparams):
        super().__init__(filename, **oparams)  # name  # Call the base class constructor

        if not oparams:
            if filename is None:
                self.velocity = params.v0
                self.gamma = params.gamma
            else:
                rows = pd.read_csv(filename)
                self.velocity = float(rows['velocity'])
                self.gamma = float(rows['gamma'])
        else:
            self.velocity = float(oparams['velocity'])
            self.gamma = float(oparams['gamma'])

        self.oparams = {'velocity': self.velocity, 'gamma': self.gamma}  # Just in case

    def calculate_effect(self, index, z, z_list, dt) -> Effect:
        """
        Effects models need a "calculate_effect" function. This function must return an Effect object
        Args:
            index (int): Index of the current enzyme
            z (Enzyme): Current enzyme
            z_list (list):
            dt (float):
        """
        # Everything 0 for now
        position = 0.0
        twist_left = 0.0
        twist_right = 0.0
        # Get neighbour enzyme
        if z.direction > 0:
            z_n = [e for e in z_list if e.position > z.position][0]  # after - On the right
        if z.direction < 0:
            z_n = [e for e in z_list if e.position < z.position][-1]  # before - On the left
        if z.direction == 0:
            raise ValueError('Error in calculating motion of RNAP. The RNAP enzyme has no direction.')

        # Check if there's one enzyme on the direction of movement. If there is one, then it will stall to avoid
        # clashing
        if z.direction > 0:  # Moving to the right
            if z_n.position - (z.position + position) <= 0:
                return position, twist_left, twist_right
        else:  # Moves to the left
            if z_n.position - (z.position + position) >= 0:
                return position, twist_left, twist_right

        # Nothing is next, so the object moves: simple uniform motion
        position, twist_left, twist_right = uniform_motion(z, dt)

        return Effect(index=index, position=position, twist_left=twist_left, twist_right=twist_right)


# ---------------------------------------------------------------------------------------------------------------------
# UTILITY FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------

# According inputs, loads the binding model, name and its params. This function is used in environment and enzyme.
# This function calls assign_effect_model
def get_effect_model(name, e_model, model_name, oparams_file, oparams):
    # If no model is given
    if e_model is None:

        # No model is given, not even a name, so there's NO effect model
        if model_name is None:
            e_model = None
            model_name = None
            oparams_file = None
            oparams = None

        # Model indicated by name
        else:
            # Loads effect model.
            # If oparams dict is given, those will be assigned to the model -> This is priority over oparams_file
            # If oparams_file is given, parameters will be read from file, in case of no oparams dict
            # If no oparams file/dict are given, default values will be used.

            # A dictionary of parameters is given so that's priority
            if isinstance(oparams, dict):
                e_model = assign_effect_model(model_name, **oparams)
            # No dictionary was given
            else:
                # If no oparams_file is given, then DEFAULT values are used.
                if oparams_file is None:
                    e_model = assign_effect_model(model_name)
                # If an oparams_file is given, then those are loaded
                else:
                    e_model = assign_effect_model(model_name, oparams_file=oparams_file)

                oparams = e_model.oparams  # To make them match

    # An actual model was given
    else:

        #  Let's check if it's actually an effect model - The model should already have the oparams
        if isinstance(e_model, EffectModel):
            #  Then, some variables are fixed.
            model_name = e_model.__class__.__name__
            oparams = e_model.oparams
            oparams_file = None

        else:
            print('Warning, effect model given is not a class for environmental ', name)
            e_model = None
            model_name = None
            oparams_file = None
            oparams = None

    return e_model, model_name, oparams_file, oparams


# Add your models into this function so it the code can recognise it
def assign_effect_model(model_name, oparams_file=None, **oparams):
    if model_name == 'RNAPUniform':
        my_model = RNAPUniform(filename=oparams_file, **oparams)
    else:
        raise ValueError('Could not recognise effect model ' + model_name)
    return my_model


def read_csv_to_dict(filename):
    """
    Reads csv file and puts it in a dictionary
    """
    return pd.read_csv(filename).to_dict()


# ---------------------------------------------------------------------------------------------------------------------
# TOPOISOMERASE ACTIVITY FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------
# Calculates the amount of coils removed by topoisomerase I
# activity. This function only depends on the supercoiling density (sigma)
# I took this function from Sam Meyer's paper (2019)
def topo1_continuum(sigma, topo, dt):
    # the function has the form of (concentration*sigmoid)*rate*dt
    a = topo.concentration * topo.k_cat * dt
    try:
        b = 1 + np.exp((sigma - topo.oparams['threshold']) / topo.oparams['width'])
        sigma_removed = a / b
    except OverflowError as oe:
        sigma_removed = 0.0
    return sigma_removed


# ----------------------------------------------------------
# Calculates the amount of coils removed by gyrase
# activity. This function only depends on the supercoiling density (sigma)
# I took this function from Sam Meyer's paper (2019)
def gyrase_continuum(sigma, gyra, dt):
    # the function has the form of (concentration*sigmoid)*rate*dt
    a = gyra.concentration * gyra.k_cat * dt
    try:
        b = 1 + np.exp(-(sigma - gyra.oparams['threshold']) / gyra.oparams['width'])
        sigma_removed = -a / b
    except OverflowError as oe:
        sigma_removed = 0.0
    return sigma_removed


# Supercoiling injection of topoisomerases. It injects according the k_cat (injected twist per second), so be careful
# because it can be both positive or negative
def topoisomerase_supercoiling_injection(topo, dt):
    position = 0.0
    # Note that k_cat is divided by two on each side because it is assumed that k_cat acts on the local region
    # (both sides)
    twist_left = 0.5 * topo.k_cat * params.w0 * dt
    twist_right = 0.5 * topo.k_cat * params.w0 * dt
    return position, twist_left, twist_right


#  def topoisomerase_lineal_supercoiling_injection(topo, dt):
#    position = 0.0


# ---------------------------------------------------------------------------------------------------------------------
# RNAP FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------

# We will use more than once this calculation, so let's store it as function
def uniform_motion(z, dt):
    # Object moves: simple uniform motion
    # position = Z.position + Z.direction * v0 * dt
    position = z.direction * v0 * dt

    # Injects twist: denatures w=gamma*v0*dt base-pairs
    twist_left = -z.direction * z.k_cat * v0 * dt
    twist_right = z.direction * z.k_cat * v0 * dt
    return position, twist_left, twist_right


# TODO: Because you have a k_cat now, this indicates how much twist RNAPs inject on each side.
# Returns new RNAP position, and the twist it injected on the left and right
def rnap_uniform_motion(z, z_list, dt):
    # Everything 0 for now
    position = 0.0
    twist_left = 0.0
    twist_right = 0.0
    # Get neighbour enzyme
    if z.direction > 0:
        z_n = [e for e in z_list if e.position > z.position][0]  # after - On the right
    if z.direction < 0:
        z_n = [e for e in z_list if e.position < z.position][-1]  # before - On the left
    if z.direction == 0:
        print('Error in calculating motion of RNAP. The RNAP enzyme has no direction.')
        sys.exit()

    # Check if there's one enzyme on the direction of movement. If there is one, then it will stall to avoid clashing
    if z.direction > 0:  # Moving to the right
        if z_n.position - (z.position + position) <= 0:
            return position, twist_left, twist_right
    else:  # Moves to the left
        if z_n.position - (z.position + position) >= 0:
            return position, twist_left, twist_right

    # Nothing is next, so the object moves: simple uniform motion
    position, twist_left, twist_right = uniform_motion(z, dt)

    return position, twist_left, twist_right


# Returns new RNAP position, and the twist it injected on the left and right - It can stall according the Geng
# model of RNAP elongation. In this model, either the RNAP moves with constant velocity or it stalls. It relies on a
# low stretching force params.f_stretching, which here we assume that all molecules interacting exert on the DNA, which
# might not be true... Additionally, if the DNA becomes hyper supercoiled, the RNAP would also stall
def rnap_torque_stall_Geng(z, z_list, dt):
    # For now, nothing happens
    position = 0.0
    twist_left = 0.0
    twist_right = 0.0
    # Get enzymes on the right and left
    z_right = [e for e in z_list if e.position > z.position][0]  # after - On the right
    z_left = [e for e in z_list if e.position < z.position][-1]  # before - On the left
    # Calculate torques and determine if the RNAP will stall
    torque_right = Marko_torque(z.superhelical)  # Torque on the right
    torque_left = Marko_torque(z_left.superhelical)  # Torque on the left
    torque = abs(torque_left - torque_right)
    if torque >= params.stall_torque:  # If torque higher than the stall torque, the RNAP stalls and doesn't move
        return position, twist_left, twist_right

    # Ok, it didn't stall, now we need the neighbours
    if z.direction > 0:
        z_n = z_right  # after - On the right
    if z.direction < 0:
        z_n = z_left  # before - On the left
    if z.direction == 0:
        print('Error in calculating motion of RNAP. The RNAP enzyme has no direction.')
        sys.exit()

    # Check if there's one enzyme on the direction of movement. If there is one, then it will stall to avoid clashing
    if z.direction > 0:  # Moving to the right
        if z_n.position - (z.position + position) <= 0:
            return position, twist_left, twist_right
    else:  # Moves to the left
        if z_n.position - (z.position + position) >= 0:
            return position, twist_left, twist_right

    # It passed all the filters and didn't stall, now, the object moves with simple uniform motion
    # position = Z.position + Z.direction * v0 * dt
    position = z.direction * v0 * dt

    # Injects twist: denatures w=gamma*v0*dt base-pairs
    twist_left = -z.direction * z.k_cat * v0 * dt
    twist_right = z.direction * z.k_cat * v0 * dt

    return position, twist_left, twist_right


# Torque calculated using Marko's elasticity model
def Marko_torque(sigma):
    if np.abs(sigma) <= np.abs(params.sigma_s):
        torque = sigma * params.cs_energy / params.w0
    elif abs(params.sigma_s) < abs(sigma) < abs(params.sigma_p):
        torque = np.sqrt(
            2 * params.p_stiffness * params.g_energy / (1 - params.p_stiffness / params.cs_energy)) / params.w0_nm
    elif abs(sigma) > abs(params.sigma_p):
        torque = sigma * params.p_stiffness / params.w0
    else:
        print('Error in Marko_torque function')
        sys.exit()
    return torque


# ---------------------------------------------------------------------------------------------------------------------
# USEFUL FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------
# This function calculates the length between two objects (proteins) considering their sizes
def calculate_length(z0, z1):
    x0 = z0.position  # positions
    x1 = z1.position
    b0 = z0.size  # size -_-
    # b1 = z1.size
    length = abs(x1 - (x0 + b0))
    # There are 4 possibilities
    # if z0.direction >= 0 and z1.direction >= 0:
    #    length = (x1 - b1) - x0
    # elif z0.direction >= 0 and z1.direction < 0:
    #    length = x1 - x0
    # elif z0.direction < 0 and z1.direction >= 0:
    #    length = (x1 - b1) - (x0 + b0)
    # elif z0.direction < 0 and z1.direction < 0:
    #    length = x1 - (x0 + b0)
    # else:
    #    print("Something went wrong in lengths")
    #    sys.exit()
    # length+=1 #1 bp needs to be added
    return length


# ----------------------------------------------------------
# This function calculates/updates the twist parameter according
# the supercoiling value of the current object Z0, and according
# the length between object Z0 and Z1.
def calculate_twist(z0, z1):
    length = calculate_length(z0, z1)  # First, I need to get the length
    sigma = z0.superhelical
    twist = sigma * w0 * length
    return twist


# ----------------------------------------------------------
# This function calculates/updates the supercoiling according
# the twist of the current object Z0, and the distance between
# Z1-Z0
def calculate_supercoiling(z0, z1):
    length = calculate_length(z0, z1)  # First, I need to get the length
    twist = z0.twist  # and twist
    if length != 0:
        sigma = twist / (w0 * length)  # and calculate the supercoiling
    else:
        sigma = 0  # I don't know if this is a solution... #But basically, this happens when a RNAP
        # binds right after one has bound
    return sigma


# ----------------------------------------------------------
# This function is equivalent to calculate_twist(), however, we use this function when
# the twist stored in the enzyme is not reliable. For example, when topoisomerases act on the DNA in the continumm
# model, we might need to convert from superhelical to twist
def calculate_twist_from_sigma(z0, z1, sigma):
    length = calculate_length(z0, z1)  # First, I need to get the length
    twist = sigma * w0 * length
    return twist


# -----------------------------------------------------------------------
# Gets the start and end positions of the fake boundaries (for circular DNA)
# In case that there is not fake boundaries, Z_N should be the last element [-1],
# in case that you have N objects including the fake boundaries, Z_N -> [N-2]
def get_start_end_c(z0, zn, nbp):
    # b_0 = z0.size
    b_n = zn.size
    x_0 = z0.position  # position of first object
    x_n = zn.position  # position of last object

    # fake position on the left
    #    position_left = 1 + x_n + b_n - nbp  # the size of the last object is considered
    position_left = x_n + b_n - nbp  # the size of the last object is considered
    # if zn.direction >= 0:  # depends on the direction
    #    position_left = 0 - (nbp - x_n)  # this is the position of the fake bit,
    # else:
    #    position_left = 0 - (nbp - (x_n + b_n))  # the size of the last object is considered

    # fake end
    position_right = nbp + x_0
    # if z0.direction >= 0:  # depends on the direction
    #    position_right = nbp + x_0 - b_0  # I think I had the sign wrong...
    # else:
    #    position_right = nbp + x_0

    return position_left, position_right

# -----------------------------------------------------------------------
