import numpy as np
from TORCphysics import params
import sys
import pandas as pd
from abc import ABC, abstractmethod

# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This module contains the mathematical functions that compose the effects model.
# Effects can describe the motion of RNAPs as well as twist injection; topoisomerase activity; overall mechanics of
# enzymes bound to DNA, etc...

# TODO: Decide which of these parameters you need
# All parameters are already in the params module, but I prefer to have them here with more simple names:
v0 = params.v0
w0 = params.w0
gamma = params.gamma


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


# Runs through each bound enzyme (in enzyme_list) and creates an effect.
# It an effects_list which contains indications on how to update the current position and how the twist on the local
# neighbouring domains are effected
def effect_model(enzyme_list, environmental_list, dt, topoisomerase_model, mechanical_model):
    # list of effects: effect = [index, position, twist_left, twist_right]
    # I use an effect list because it's easier because there are multiple changes in the local twists
    effect_list = []
    for i, enzyme in enumerate(enzyme_list):

        if enzyme.enzyme_type == 'EXT':  # We can speed up things a bit by ignoring the fake boundaries
            continue

        # TODO: add effect for other type of enzymes, like NAPs? Maybe their effect is to do nothing. And definitely
        #  an effect for topos
        # Administer the effect model to use.
        # -------------------------------------------------------------------------------------------------------------
        # From these models, they can update the position of the current enzyme and affect the local twist on the right
        # and left
        if enzyme.enzyme_type == 'RNAP':  # For now, only RNAPs have an effect
            # TODO: in the future, according the input we may choose between different motion models, maybe one with
            #  torques and not uniform motion
            if mechanical_model == 'uniform':
                # Calculates change in position and the twist that it injected on the left and right
                position, twist_left, twist_right = rnap_uniform_motion(enzyme, enzyme_list, dt)
            elif mechanical_model == 'torque_stall_Geng':
                position, twist_left, twist_right = rnap_torque_stall_Geng(enzyme, enzyme_list, dt)
            else:
                print('Sorry, cannot recognise mechanistic model')
                sys.exit()
            # size = abs(enzyme.site.start - enzyme.site.end + 1)
            # output_environment = Environment(e_type='mRNA', name=enzyme.site.name, site_list=[], concentration=1,
            #                                 k_on=0, k_off=0, k_cat=0, size=size)
        #            output_enzyme = Enzyme(e_type='mRNA', name=enzyme.site.name, site=None, position=None, size=size,
        #                                   twist=0, superhelical=0)
        elif enzyme.enzyme_type == 'topo':
            # topo = [environment for environment in environmental_list
            #         if environment.name == enzyme.name][0]  # Can select the model from here?
            position, twist_left, twist_right = topoisomerase_supercoiling_injection(enzyme, dt)
            # TODO: It would probably be better if we have sub-models for effects?
            #  So we would have a binding model and a effect model.
        #           if topoisomerase_model == 'stochastic':
        #                position, twist_left, twist_right = topoisomerase_supercoiling_injection(enzyme, dt)
        #            if topoisomerase_model == 'random_lineal':
        #                position, twist_left, twist_right = topoisomerase_lineal_supercoiling_injection(enzyme, dt)
        else:
            continue

        # Now create the effect taken place at the enzyme i
        effect_list.append(Effect(index=i, position=position, twist_left=twist_left, twist_right=twist_right))

    # Topoisomerase continuum model - If we are using a continuum model, then we need to add the topos effects
    # --------------------------------------------------------------
    if topoisomerase_model == 'continuum':
        # Gets list of topoisomerase enzymes in the environment
        topo_list = [environment for environment in environmental_list
                     if environment.enzyme_type == 'topo' or environment.enzyme_type == 'topoisomerase']
        position = 0  # topoisomerases cannot change enzymes positions and only affect each local site (twist_right)
        twist_left = 0
        for topo in topo_list:
            for i, enzyme in enumerate(enzyme_list):
                # We can speed up things a bit by ignoring the fake boundaries
                if enzyme.name == 'EXT_L' and len(enzyme_list) > 2:
                    continue
                elif enzyme.name == 'EXT_R':
                    continue

                if topo.name == 'topoI':
                    sigma = topo1_continuum(enzyme.superhelical, topo, dt)
                    twist_right = calculate_twist_from_sigma(enzyme, enzyme_list[i + 1], sigma)
                elif topo.name == 'gyrase':
                    sigma = gyrase_continuum(enzyme.superhelical, topo, dt)
                    twist_right = calculate_twist_from_sigma(enzyme, enzyme_list[i + 1], sigma)
                else:
                    twist_right = 0  # If it doesn't recognize the topo, then don't do anything

                # Now create the effect taken place at the enzyme i
                effect_list.append(Effect(index=i, position=position, twist_left=twist_left, twist_right=twist_right))

    return effect_list


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

# ---------------------------------------------------------------------------------------------------------------------
# EFFECT MODELS
# ---------------------------------------------------------------------------------------------------------------------
# Add your models into this function so it the code can recognise it
def assign_effect_model(model_name, oparams_file, enzyme_name):
    if model_name == 'RNAPUniform':
        my_model = RNAPUniform()
    else:
        sys.exit()
    return my_model


class EffectModel(ABC):
    """
     The EffectModel class for defining effect models.
     If you need a new one, define it here.
    """

    #   def __init__(self, name):
    #        """
    #        Args:
    #            name (str): EffectModels should have a name as a minimum.
    #        """
    #        self.name = self.__class__.__name__

    def __init__(self):
        pass

    @abstractmethod
    def calculate_effect(self) -> Effect:
        """
        Effects models need a "calculate_effect" function. This function must return an Effect object
        """
        pass


class RNAPUniform(EffectModel):
    """
     Effect model for the RNAP uniform motion. RNAP moves with constant speed and injects +-supercoils on each site
    """

    # def __init__(self, name, filename):
    def __init__(self, filename):

#            """
#        Args:
#            name (str): Name
#        Optional:
#            filename (str): Path to parametrization of the model. If not specified, default values will be used.
#        """
        super().__init__()  # Call the base class constructor

        if filename is None or filename == 'none':
            self.velocity = params.v0
            self.gamma = params.gamma
        else:
            oparams = read_csv_to_dict(filename)
            self.velocity = oparams['velocity']
            self.gamma = oparams['gamma']

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
            print('Error in calculating motion of RNAP. The RNAP enzyme has no direction.')
            sys.exit()

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


def read_csv_to_dict(filename):
    """
    Reads csv file and puts it in a dictionary
    """
    return pd.read_csv(filename).to_dict()
