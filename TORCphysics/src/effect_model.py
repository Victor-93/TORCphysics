import numpy as np
from TORCphysics import params
import sys

# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This module contains the mathematical functions that compose the effects model.
# Effects can describe the motion of RNAPs as well as twist injection; topoisomerase activity; overall mechanics of
# enzymes bound to DNA, etc...

# All parameters are already in the params module, but I prefer to have them here with more simple names:
v0 = params.v0
w0 = params.w0
gamma = params.gamma
topo_w = params.topo_w
topo_t = params.topo_t
gyra_w = params.gyra_w
gyra_t = params.gyra_t


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
                position, twist_left, twist_right = rnap_uniform_motion(enzyme, dt)
            else:
                print('Sorry, cannot recognise mechanistic model')
                sys.exit()
            # size = abs(enzyme.site.start - enzyme.site.end + 1)
            # output_environment = Environment(e_type='mRNA', name=enzyme.site.name, site_list=[], concentration=1,
            #                                 k_on=0, k_off=0, k_cat=0, size=size)
        #            output_enzyme = Enzyme(e_type='mRNA', name=enzyme.site.name, site=None, position=None, size=size,
        #                                   twist=0, superhelical=0)
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
                    sigma = topo1_continuum(enzyme.superhelical, topo.concentration, topo.k_cat, dt)
                    twist_right = calculate_twist_from_sigma( enzyme, enzyme_list[i+1], sigma )
                elif topo.name == 'gyrase':
                    sigma = gyrase_continuum(enzyme.superhelical, topo.concentration, topo.k_cat, dt)
                    twist_right = calculate_twist_from_sigma( enzyme, enzyme_list[i+1], sigma )
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
def topo1_continuum(sigma, topo_c, topo_k, dt):
    # the function has the form of (concentration*sigmoid)*rate*dt
    a = topo_c * topo_k * dt
    try:
        b = 1 + np.exp((sigma - topo_t) / topo_w)
        sigma_removed = a / b
    except OverflowError as oe:
        sigma_removed = 0.0
    return sigma_removed


# ----------------------------------------------------------
# Calculates the amount of coils removed by gyrase
# activity. This function only depends on the supercoiling density (sigma)
# I took this function from Sam Meyer's paper (2019)
def gyrase_continuum(sigma, gyra_c, gyra_k, dt):
    # the function has the form of (concentration*sigmoid)*rate*dt
    a = gyra_c * gyra_k * dt
    try:
        b = 1 + np.exp(-(sigma - gyra_t) / gyra_w)
        sigma_removed = -a / b
    except OverflowError as oe:
        sigma_removed = 0.0
    return sigma_removed


# Returns new RNAP position, and the twist it injected on the left and right
def rnap_uniform_motion(z, dt):
    # Object moves: simple uniform motion
    # position = Z.position + Z.direction * v0 * dt
    position = z.direction * v0 * dt
    # Injects twist: denatures w=gamma*v0*dt base-pairs
    twist_left = -z.direction * gamma * v0 * dt
    twist_right = z.direction * gamma * v0 * dt
    return position, twist_left, twist_right


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
