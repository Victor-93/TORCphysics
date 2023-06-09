import sys

import numpy as np
from TORCphysics import params, Enzyme

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
topo_w = params.topo_w
topo_t = params.topo_t
gyra_w = params.gyra_w
gyra_t = params.gyra_t


# TODO: needs a unbinding model selector.
# TODO: I need to find a way to set a method so it is easy to define binding models
# TODO: We still need to check if the new enzymes are not overlapping. If more than 1 enzyme
#  passed the probability test and are binding the same overlapping region, we need to flip a coin to
#  check which one will be binding

# ---------------------------------------------------------------------------------------------------------------------
# BINDING/UNBINDING MODEL SELECTION FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------
# If you create a new model, modify these functions so the program can process them, and select the appropriate
# binding and unbinding models

# This function administrates the binding model to use
# ---------------------------------------------------------------------------------------------------------------------
def select_binding_model(site, environment, site_superhelical, dt):
    have_model = True  # Tells the function that called if the model exists
    rate = np.zeros_like(site_superhelical)
    binding_probability = 0.0

    # Simple poisson process (constant binding)
    if site.site_model == 'poisson' or site.site_model == 'Poisson':
        rate = site.k_min * np.ones_like(site_superhelical)
        binding_probability = P_binding_Poisson(site.k_min * np.ones_like(site_superhelical), dt)
    # MODELS - This models include all enzymes?:
    # Sam's Meyer model
    elif site.site_model == 'sam' or site.site_model == 'Sam':
        rate = promoter_curve_Meyer(site.k_min, site_superhelical)
        binding_probability = P_binding_Nonh_Poisson(rate, dt)
    # Max-min model according oparams measured with SIDD
    elif site.site_model == 'maxmin' or site.site_model == 'Maxmin':
        rate = promoter_curve_opening_E_maxmin(site.k_min, site.k_max, site_superhelical, *site.oparams)
        binding_probability = P_binding_Nonh_Poisson(rate, dt)
    # Inverted max-min model (where it is positive supercoiling sensitive)
    elif site.site_model == 'maxmin_I' or site.site_model == 'Maxmin_I':
        rate = promoter_curve_opening_E_maxmin_I(site.k_min, site.k_max, site_superhelical, *site.oparams)
        binding_probability = P_binding_Nonh_Poisson(rate, dt)
    # Similar to max-min but with the effective energy
    elif site.site_model == 'effE' or site.site_model == 'EffE':
        rate = promoter_curve_opening_E(site.k_min, site_superhelical, sigma0=0, *site.oparams)
        binding_probability = P_binding_Nonh_Poisson(rate, dt)
    elif site.site_model == 'none' or site.site_model == 'None' or site.site_model is None:
        have_model = False
    elif site.site_model == 'stochastic_topoI':
        rate = topoI_binding(environment, site_superhelical)
        binding_probability = P_binding_Nonh_Poisson(rate, dt)
    elif site.site_model == 'stochastic_gyrase':
        rate = gyrase_binding(environment, site_superhelical)
        binding_probability = P_binding_Nonh_Poisson(rate, dt)
    else:  # If there's no model, there's no binding
        have_model = False
    return rate, binding_probability, have_model


# This function is in charge of administrating the unbinding models to use
# ---------------------------------------------------------------------------------------------------------------------
def select_unbinding_model(enzyme, dt, rng):
    unbind = False
    if enzyme.enzyme_type == 'RNAP':
        have_model = True
        unbind = RNAP_unbinding_model(enzyme)
    # TODO: some NAPs might not follow a Poisson_unbinding_model
    elif enzyme.enzyme_type == 'topo' or enzyme.enzyme_type == 'NAP':
        have_model = True
        unbind = Poisson_unbinding_model(enzyme, dt, rng)
    else:
        # TODO: When you write the warnings, add this one. And do something similar for the effects model
        # print('Warning, we do not have the unbinding model for your enzyme type:', enzyme.enzyme_type)
        have_model = False
    return unbind, have_model


# ---------------------------------------------------------------------------------------------------------------------
# BINDING/UNBINDING OVERALL MODELS
# ---------------------------------------------------------------------------------------------------------------------
# These functions do not perform a particular binding model, but are the ones that implement the binding/unbinding
# processes.


# Goes through enzymes in the environment (environmental_list) and search for the available sites that it recognizes.
# If the site is available, then, according site model (and in the future also environmental model) calculate the
# binding probability. It then returns a list of new_enzymes that will bind the DNA
# rng - is a numpy random generator
def binding_model(enzyme_list, environmental_list, dt, rng):
    new_enzymes = []  # here I will include the new enzymes

    # Go through environment
    for i, environment in enumerate(environmental_list):

        # For now, only RNAPs/genes
        # if environment.enzyme_type != 'RNAP':
        #    continue

        # If we ran out of the enzyme in the environment, then there's nothing to do
        if environment.concentration <= 0.0:
            continue

        # Go through sites
        for j, site in enumerate(environment.site_list):
            # We will do the binding process in this order:
            # Check site and model to use.
            # And calculate binding probability and if it'll bind
            # Check if binding site is available
            # If there are multiple enzymes that want to bind but their ranges overlap, we must choose

            # For now, only genes!
            # -----------------------------------------------------------
            # if site.site_type != 'gene':
            #    continue
            if '_global' in site.name:  # We don't actually model it for globals
                continue

            # Get superhelical density at site
            enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= site.start][-1]
            site_superhelical = enzyme_before.superhelical

            # According model, calculate the binding probability
            # -----------------------------------------------------------
            # We need to figure out 1 or 2 models.
            # A model for the rate in case it is modulated by supercoiling
            # And a model for calculating the binding probability.

            rate, binding_probability, have_model = select_binding_model(site, environment, site_superhelical, dt)
            if not have_model:  # If I don't have a model, then we skip
                continue

            # Decide if the enzyme will bind
            # -------------------------------------------------------------
            urandom = rng.uniform()  # we need a random number
            # print(environment.name, urandom, binding_probability)

            if urandom <= binding_probability:  # and decide

                # Check if site is available  - it is actually faster to first calculate the probability, so I move
                # it here.
                # -------------------------------------------------------------
                site_available = check_site_availability(site, enzyme_list, environment.size)
                if not site_available:
                    continue

                # Add enzyme
                # --------------------------------------------------------
                # We first need to figure out the position, twist and superhelical (these last two will be sorted in
                # the circuit module

                # position:
                if site.direction > 0:
                    position = site.start - environment.size
                elif site.direction <= 0:
                    position = site.start
                else:
                    print("Error in adding new enzyme to list of new enzymes")
                    sys.exit()

                # Create enzyme, and note that it is missing twist and the superhelical density.
                # I think it's better to fix it in the circuit module
                enzyme = Enzyme(e_type=environment.enzyme_type, name=environment.name, site=site, position=position,
                                size=environment.size, k_cat=environment.k_cat, k_off=environment.k_off,
                                twist=0.0, superhelical=0.0)

                new_enzymes.append(enzyme)

    # TODO: check_binding_conflicts() needs testing
    new_enzymes = check_binding_conflicts(new_enzymes, rng)

    return new_enzymes


# Goes through the enzymes in enzymes list and according to their unbinding condition unbind them.
# Returns a list of enzyme indexes that will unbind, the enzyme that unbinds
# ---------------------------------------------------------------------------------------------------------------------
def unbinding_model(enzymes_list, dt, rng):
    drop_list_index = []  # This list will have the indices of the enzymes that will unbind, and the enzyme
    drop_list_enzyme = []  # And a list with the enzymes
    for i, enzyme in enumerate(enzymes_list):

        if enzyme.enzyme_type == 'EXT':  # The fake boundaries can't unbind
            continue

        # According enzyme_type, apply unbinding condition
        # ------------------------------------------------------------------
        unbind, have_model = select_unbinding_model(enzyme, dt, rng)
        if not have_model:  # If I don't have a model, then we skip
            continue

        # Now add it to the drop_list if the enzyme will unbind
        # ------------------------------------------------------------------
        if unbind:
            drop_list_index.append(i)
            drop_list_enzyme.append(enzyme)

    return drop_list_index, drop_list_enzyme


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
def Poisson_unbinding_model(enzyme, dt, rng):
    unbinding_probability = P_binding_Poisson(enzyme.k_off, dt)
    urandom = rng.uniform()  # we need a random number
    if urandom <= unbinding_probability:  # and decide
        unbind = True
    else:
        unbind = False
    return unbind


# ---------------------------------------------------------------------------------------------------------------------
# BINDING MODELS
# ---------------------------------------------------------------------------------------------------------------------
# These functions are for the particular binding model.
# If you need a new one, define it here.

# ----------------------------------------------------------
# This equation calculates the probability of binding according
# the Poisson process
def P_binding_Poisson(rate, dt):
    rdt = rate * dt  # it is what is in the exponent (is that how you call it?)
    probability = rdt * np.exp(-rdt)

    return probability


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
    b = 1 + np.exp((sigma - topo.oparams['threshold']) / topo.oparams['width']) # [0] because the dictionary
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
                del enzyme_list[i-s-1]
                s += 1
                # checked_enzyme_list.append(enzyme_before)
            else:               # And if >0.5, then we don't add the enzyme before (we lose it).
                del enzyme_list[i-s]
                s += 1
                # continue
        # else:
            # checked_enzyme_list.append(enzyme_before)  # If nothing overlaps, then nothing happens


    return enzyme_list  #checked_enzyme_list
