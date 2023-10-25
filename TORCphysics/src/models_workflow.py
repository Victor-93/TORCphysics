import numpy as np
from TORCphysics import utils, Enzyme
import sys


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
    """
    This function implements the binding workflow, where a list of newly bound enzymes are returned.
    Basically, it goes through the environmentals in environmental_list and the sites they recognise to determine
    if they will ind each of them. Some binding models depend on local variables such as the superhelical density,
    which are stored in the Enzymes in the enzyme_list.

    Inputs
    ----------
    enzyme_list : list
        This is a list of Enzymes that are currently bound to the DNA.
    environmental_list : list
        This is a list of Environmentals.
    dt : float
        Simulation timestep in seconds (s).
    rng : numpy object
        Random number generator. This object needs to be previously defined.

    Returns
    ----------
    new_enzymes : list
        This is a list that contains the new Enzymes to be added to the current enzyme_list.
    """

    new_enzymes = []  # here I will include the new enzymes

    # We will do the binding process in this order:
    # Go through the sites each environmental recognise.
    # For each of those sites, check the model to use.
    # Then, calculate the binding probability and if the enzyme will bind, check if the binding site is available
    # If there are multiple enzymes that want to bind but their ranges overlap, we must choose which will bind.

    # Go through environment
    for i, environment in enumerate(environmental_list):

        # If we ran out of the enzyme in the environment, then there's nothing to do
        if environment.concentration <= 0.0:
            continue

        # Go through sites
        for j, site in enumerate(environment.site_list):

            # Pass if site does not have binding model
            if site.binding_model is None:
                continue

            if site.global_site:  # We don't actually model it for globals (they can't be bind)
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
                site_available = utils.check_site_availability(site=site, enzyme_list=enzyme_list,
                                                               environmental=environment)
                if not site_available:
                    continue

                # Add enzyme
                # --------------------------------------------------------
                # We still need to figure out the position, twist and superhelical, but these last two will be sorted in
                # the circuit module

                # TODO: Make sure that this applies to any sites and enzymes binding.
                #  It means that enzymes always bind on the behind the start site. This makes sense for genes.
                #  I think this is correct! Draw it, think about it, and apply it everywhere.
                #  And think about how TF would work, or repressors, they should block the sites, but TF should
                #  increase the binding probability.
                # TODO: Think about this,
                # position:
                if site.direction > 0:
                    position = site.start - environment.size
                elif site.direction <= 0:
                    position = site.start
                else:
                    print("Error in adding new enzyme to list of new enzymes")
                    sys.exit()

                # Create enzyme, and note that it is missing twist and the superhelical density.
                # Those will be added in the circuit module
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
    """
    This function implements the effect workflow, where a list of Effects that affect the local variables of DNA,
    are returned.
    Basically, it goes through the list of bound Enzymes stored in enzyme_list and implement the EffectModel of each
    of them.
    Some Environmentals have a continuum representation of their activity, and those can affect every local region in
    the DNA defined by the bound enzymes.

    Inputs
    ----------
    enzyme_list : list
        This is a list of Enzymes that are currently bound to the DNA.
    environmental_list : list
        This is a list of Environmentals.
    dt : float
        Simulation timestep in seconds (s).

    Returns
    ----------
    effect : list
        This is a list that contains the Effects that will be applied to each local region on the DNA.
    """

    # list of effects: effect = [index, position, twist_left, twist_right]
    # I use an effect list because it's easier because there are multiple changes in the local twists
    effect_list = []

    # The plan is to go through each enzyme in enzyme_list, and apply their effect_model
    for i, enzyme in enumerate(enzyme_list):

        if enzyme.enzyme_type == 'EXT':  # We can speed up things a bit by ignoring the fake boundaries
            continue

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

    # Topoisomerase continuum model - If we are using a continuum model, then we need to add the topos effects
    # --------------------------------------------------------------
    # The plan is to go through the list of environmentals in environmental_list, abd if they have a continuum
    # effect_model, then go through the list of enzymes and apply the effect on each of the local regions.
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

    return effect_list


# Unbinding Workflow
# ---------------------------------------------------------------------------------------------------------------------

# Goes through the enzymes in enzymes list and according to their unbinding condition unbind them.
# Returns a list of enzyme indexes that will unbind, the enzyme that unbinds
# ---------------------------------------------------------------------------------------------------------------------
# TODO: Think if substances to the environment will be realised with these reactions, e.g., maybe mRNA will be realised?
def unbinding_workflow(enzymes_list, dt, rng):
    """
    This function implements the unbinding workflow, where a list of unbinding enzymes is returned.
    Basically, it goes through the enzymes in enzymes_list and according to their unbinding probability calculated with
    their unbinding_model, they will unbind the DNA.

    Inputs
    ----------
    enzyme_list : list
        This is a list of Enzymes that are currently bound to the DNA.
    dt : float
        Simulation timestep in seconds (s).
    rng : numpy object
        Random number generator. This object needs to be previously defined.

    Returns
    ----------
    drop_list_index : list
        This is a list with the indices of the enzymes (in enzyme_list) that will unbind.
    drop_list_enzyme : list
        And a list of the same enzymes that will unbind the DNA.
    """

    drop_list_index = []  # This list will have the indices of the enzymes that will unbind, and the enzyme
    drop_list_enzyme = []  # And a list with the enzymes

    # Go through each enzyme and determine if they will unbind.
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

    return drop_list_index, drop_list_enzyme


# ---------------------------------------------------------------------------------------------------------------------
# HELPFUL FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------
# TODO: Maybe in the future, when including TF, this function could be included in the binding model.
#  But for the mean time, just leave it here. Maybe TFs don't work like I have on my mind, which is that they
#  provide an adhesive energy and in consequence, increase the rate of binding of RNAPs, and when RNAPs bind, the TF
#  unbind? Maybe this isn't the case.

# site = site to be bound
# environmental = enzyme binding the site
# enzyme_list = list of enzymes currently bound to the DNA
def check_site_availability3(site, environmental, enzyme_list):

    # Enzyme before the start site
    e_b = [enzyme for enzyme in enzyme_list if enzyme.position <= site.start][-1]
    # Enzyme after the start site
    e_a = [enzyme for enzyme in enzyme_list if enzyme.position >= site.start][0]

    # Let's build ranges
    # ----------------------------------------------------
    # Range of binding enzyme. This the region that the environmental is trying to bind in the DNA
    bind_r1 = site.start - (environmental.size - environmental.effective_size)*.05
    bind_r2 = site.start + (environmental.size + environmental.effective_size)*.05
    # Range of DNA that occupies the enzyme before
    r_b1 = e_b.position - (e_b.size - e_b.effective_size)*0.5
    r_b2 = e_b.position + (e_b.size + e_b.effective_size)*0.5
    # Range of DNA that occupies the enzyme after
    r_a1 = e_a.position - (e_a.size - e_a.effective_size)*0.5
    r_a2 = e_a.position + (e_a.size + e_a.effective_size)*0.5

    # Check if ranges overlap
    # ----------------------------------------------------

    # Let's check for the enzyme on the left (r_b)
    overlap = max(r_b1, bind_r1) <= min(r_b2, bind_r2)

    if overlap:
        available = False
        return available

    # Let's check for the enzyme on the right (r_a)
    overlap = max(bind_r1, r_a1) <= min(bind_r2, r_a2)

    if overlap:
        available = False
        return available

    # TODO: I did something, but does the direction matter? yes! Check how!
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

    # TODO: The function does not work!
    # If any of them intersect
    if (set(range_before) & set(my_range)) or (set(range_after) & set(my_range)):
        available = False
    # there is an intersection
    else:
        available = True

    return available



# This function checks if a site is not blocked by other enzymes
def check_site_availability2(site, enzyme_list, size):
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

    # TODO: The function does not work!
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
# binding3
# TODO: test this function - design a experiment in which you kind of know what outcome you should get.
# TODO: pass this enzye to utils
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
