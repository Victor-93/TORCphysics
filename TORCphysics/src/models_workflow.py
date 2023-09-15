import numpy as np
from TORCphysics import params
import sys
import pandas as pd
from abc import ABC, abstractmethod


# Runs through each bound enzyme (in enzyme_list) and creates an effect.
# It an effects_list which contains indications on how to update the current position and how the twist on the local
# neighbouring domains are effected
def effect_workflow(enzyme_list, environmental_list, dt, topoisomerase_model, mechanical_model):
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
    elif 'poisson_lineal' in site.site_model:
        rate = environment.k_on * np.ones_like(site_superhelical)
        binding_probability = P_binding_Poisson(rate, dt)
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
