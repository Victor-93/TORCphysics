import numpy as np
from TORCphysics import params, Enzyme
import sys
import pandas as pd
from abc import ABC, abstractmethod


# TODO: Properly sort, comment and document your functions.

# ---------------------------------------------------------------------------------------------------------------------
# GENERAL WORKFLOWS
# ---------------------------------------------------------------------------------------------------------------------
# These functions do not perform a particular binding model, but are the ones that implement the binding/unbinding
# processes.


# Binding Workflow
# ---------------------------------------------------------------------------------------------------------------------

# Goes through enzymes in the environment (environmental_list) and search for the available sites that it recognizes.
# If the site is available, then, according site binding model calculate the binding probability.
# It then returns a list of new_enzymes that will bind the DNA
# rng - is a numpy random generator
def binding_workflow(enzyme_list, environmental_list, dt, rng):
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

            # Pass if site does not have binding model
            if site.binding_model is None:
                continue
            # We will do the binding process in this order:
            # Check site and model to use.
            # And calculate binding probability and if it'll bind
            # Check if binding site is available
            # If there are multiple enzymes that want to bind but their ranges overlap, we must choose

            # For now, only genes!
            # -----------------------------------------------------------
            # if site.site_type != 'gene':
            #    continue
            # TODO: Use site_global
            if '_global' in site.name:  # We don't actually model it for globals
                continue

            # Get superhelical density at site
            enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= site.start][-1]
            site_superhelical = enzyme_before.superhelical

            # According model, calculate the binding probability
            # -----------------------------------------------------------
            binding_probability = site.binding_model.binding_probability(environmental=environment,
                                                                         superhelical=site_superhelical, dt=dt)

            # Decide if the enzyme will bind
            # -------------------------------------------------------------
            urandom = rng.uniform()  # we need a random number

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

                # TODO: Make sure that this applies to any sites and enzymes binding.
                #  It means that enzymes always bind on the left from the start site.
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

                # TODO: Check the unbinding model and effect model. Some unbinding models depend on the site parameters,
                #  for example, in genes, their unbinding model depends on the site. But the site is already provided,
                #  so maybe is better that the unbinding model does this on its own.
                enzyme = Enzyme(e_type=environment.enzyme_type, name=environment.name, site=site, position=position,
                                size=environment.size, effective_size=environment.effective_size, twist=0.0,
                                superhelical=0.0, unbinding_model=environment.unbinding_model)
                new_enzymes.append(enzyme)

    # TODO: check_binding_conflicts() needs testing
    new_enzymes = check_binding_conflicts(new_enzymes, rng)

    # TODO: In the future, It can also return changes in the environment. Maybe the binding of an enzyme changes the
    #  concentration? Or not...?
    return new_enzymes


# Effect Workflow
# ---------------------------------------------------------------------------------------------------------------------

# Runs through each bound enzyme (in enzyme_list) and creates an effect.
# It an effects_list which contains indications on how to update the current position and how the twist on the local
# neighbouring domains are effected
def effect_workflow(enzyme_list, environmental_list, dt):
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

        # Check if bound enzyme has effect model
        if enzyme.effect_model is None:
            continue
        # Calculate effect and add it to the list
        effect_i = enzyme.effect_model.calculate_effect(index=i, z=enzyme, z_list=enzyme_list, dt=dt)
        effect_list.append(effect_i)

        # if enzyme.enzyme_type == 'RNAP':  # For now, only RNAPs have an effect
        #    # TODO: in the future, according the input we may choose between different motion models, maybe one with
        #    #  torques and not uniform motion
        #    if mechanical_model == 'uniform':
        #        # Calculates change in position and the twist that it injected on the left and right
        #        position, twist_left, twist_right = rnap_uniform_motion(enzyme, enzyme_list, dt)
        #    elif mechanical_model == 'torque_stall_Geng':
        #        position, twist_left, twist_right = rnap_torque_stall_Geng(enzyme, enzyme_list, dt)
        #    else:
        #        print('Sorry, cannot recognise mechanistic model')
        #        sys.exit()
        #    # size = abs(enzyme.site.start - enzyme.site.end + 1)
        #    # output_environment = Environment(e_type='mRNA', name=enzyme.site.name, site_list=[], concentration=1,
        #    #                                 k_on=0, k_off=0, k_cat=0, size=size)
        #            output_enzyme = Enzyme(e_type='mRNA', name=enzyme.site.name, site=None, position=None, size=size,
        #                                   twist=0, superhelical=0)
        # elif enzyme.enzyme_type == 'topo':
        #    # topo = [environment for environment in environmental_list
        #    #         if environment.name == enzyme.name][0]  # Can select the model from here?
        #    position, twist_left, twist_right = topoisomerase_supercoiling_injection(enzyme, dt)
        # else:
        #    continue

        # Now create the effect taken place at the enzyme i
        # effect_list.append(Effect(index=i, position=position, twist_left=twist_left, twist_right=twist_right))

    # Topoisomerase continuum model - If we are using a continuum model, then we need to add the topos effects
    # --------------------------------------------------------------
    for environmental in environmental_list:
        if environmental.effect_model is None:
            continue
        if environmental.effect_model.continuum:
            for i, enzyme in enumerate(enzyme_list):

                # We can speed up things a bit by ignoring the fake boundaries
                if enzyme.name == 'EXT_L' and len(enzyme_list) > 2:
                    continue
                elif enzyme.name == 'EXT_R':
                    continue

                # Calculate effect and add it to the list
                effect_i = environmental.effect_model.calculate_effect(concentration=environmental.concentration,
                                                                       index=i, z=enzyme, z_list=enzyme_list, dt=dt)
                effect_list.append(effect_i)

#    if topoisomerase_model == 'continuum':
#        # Gets list of topoisomerase enzymes in the environment
#        topo_list = [environment for environment in environmental_list
#                     if environment.enzyme_type == 'topo' or environment.enzyme_type == 'topoisomerase']
#        position = 0  # topoisomerases cannot change enzymes positions and only affect each local site (twist_right)
#        twist_left = 0
#        for topo in topo_list:
#            for i, enzyme in enumerate(enzyme_list):
#                # We can speed up things a bit by ignoring the fake boundaries
#                if enzyme.name == 'EXT_L' and len(enzyme_list) > 2:
#                    continue
#                elif enzyme.name == 'EXT_R':
#                    continue

#                if topo.name == 'topoI':
#                    sigma = topo1_continuum(enzyme.superhelical, topo, dt)
#                    twist_right = calculate_twist_from_sigma(enzyme, enzyme_list[i + 1], sigma)
#                elif topo.name == 'gyrase':
#                    sigma = gyrase_continuum(enzyme.superhelical, topo, dt)
#                    twist_right = calculate_twist_from_sigma(enzyme, enzyme_list[i + 1], sigma)
#                else:
#                    twist_right = 0  # If it doesn't recognize the topo, then don't do anything

#                # Now create the effect taken place at the enzyme i
#                effect_list.append(Effect(index=i, position=position, twist_left=twist_left, twist_right=twist_right))

    return effect_list


# Unbinding Workflow
# ---------------------------------------------------------------------------------------------------------------------

# Goes through the enzymes in enzymes list and according to their unbinding condition unbind them.
# Returns a list of enzyme indexes that will unbind, the enzyme that unbinds
# ---------------------------------------------------------------------------------------------------------------------
def unbinding_workflow(enzymes_list, dt, rng):
    drop_list_index = []  # This list will have the indices of the enzymes that will unbind, and the enzyme
    drop_list_enzyme = []  # And a list with the enzymes
    for i, enzyme in enumerate(enzymes_list):

        if enzyme.enzyme_type == 'EXT':  # The fake boundaries can't unbind
            continue

        if enzyme.unbinding_model is None:  # Skip if enzyme doesn't have unbinding model (cannot unbind)
            continue

        # According enzyme_type, apply unbinding condition
        # ------------------------------------------------------------------
        unbinding_probability = enzyme.unbinding_model.unbinding_probability(enzyme, dt)

        urandom = rng.uniform()  # we need a random number

        if urandom <= unbinding_probability:  # and decide if it'll unbind
            drop_list_index.append(i)
            drop_list_enzyme.append(enzyme)

        #  unbind, have_model = select_unbinding_model(enzyme, dt, rng)
        #  if not have_model:  # If I don't have a model, then we skip
         #   continue

        # Now add it to the drop_list if the enzyme will unbind
        # ------------------------------------------------------------------
        #  if unbind:
        #    drop_list_index.append(i)
        #    drop_list_enzyme.append(enzyme)

    return drop_list_index, drop_list_enzyme



# ---------------------------------------------------------------------------------------------------------------------
# HELPFUL FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------
# If you create a new model, modify these functions so the program can process them, and select the appropriate
# binding and unbinding models

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


# ---------------------------------------------------------------------------------------------------------------------
# OBSOLETE FUNCTIONS (Erase them when you can)
# ---------------------------------------------------------------------------------------------------------------------

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


