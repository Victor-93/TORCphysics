import numpy as np
from TORCphysics import params, Enzyme, Environment

# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This module contains the mathematical functions that compose the statistical part of my model
# It comprasis the necessary equations required for simulating the stochastic binding, and the
# topoisomerases/gyrases activities

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


#class Add_enzyme:
#    # For handling the data in the data frames, it is more useful to create a class that contains information
#    # regarding the twist and superhelical density at the site prior binding.#

#    # So, the twist and superhelical density are the local values at the site prior binding.
#    def __init__(self, will_bind, enzyme, twist, superhelical):
#        # I'll save the input filenames just in case
#        self.will_bind = will_bind
#        self.enzyme = enzyme
#        self.twist = twist
#        self.superhelical = superhelical

# ---------------------------------------------------------------------------------------------------------------------
# BINDING FUNCTIONS
# ---------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------

# This function is in charge of administrating the binding model.
# According to the site's type and properties, it determines which binding model and probabilistic model to use
# and then uses those models to determine if enzymes in the environment would actually bind the DNA

# probabilistic model we will use for binding events.
# model = poisson or not (non homogeneus Poisson)
# genome = the dataframe with all the information in the piece of "genome"
# sigma  = supercoiling density that the binding site
# f_ini  = initiation function

# Goes through the enzymes in enzymes list and according to their unbinding condition unbind them.
# Returns a list of enzyme indexes that will unbind, the enzyme that unbinds
def unbinding_model(enzymes_list):
    drop_list_index = []  # This list will have the indices of the enzymes that will unbind, and the enzyme
    drop_list_enzyme = []  # And a list with the enzymes
    for i, enzyme in enumerate(enzymes_list):

        if enzyme.enzyme_type == 'EXT':  # The fake boundaries can't unbind
            continue

        unbind = False  # True if it unbinds

        # According enzyme_type, apply unbinding condition
        # ------------------------------------------------------------------
        # TODO: add unbinding condition for more enzyme types!
        # Right now, only RNAPs can unbind
        if enzyme.enzyme_type == 'RNAP':
            # condition for transcription in >>>>>>>>>>>>> right direction or
            # condition for transcription in <<<<<<<<<<<<< left  direction
            if (enzyme.direction == 1 and enzyme.end - enzyme.position <= 0) or \
                    (enzyme.direction == -1 and enzyme.end - enzyme.position >= 0):
                unbind = True
        else:
            continue

        # Now add the to the drop_list if the enzyme will unbind
        # ------------------------------------------------------------------
        if unbind:
            drop_list_index.append(i)
            drop_list_enzyme.append(enzyme)

    return drop_list_index, drop_list_enzyme


# Goes through enzymes in the environment (environmental_list) and search for the available sites that it recognizes.
# If the site is available, then, according site model (and in the future also environmental model) calculate the
# binding probability. It then returns a list of new_enzymes that will bind the DNA
def binding_model(enzyme_list, environmental_list, dt):
    new_enzymes = []  # here I will include the new enzymes

    # Go through environment
    for i, environment in enumerate(environmental_list):

        # For now, only RNAPs/genes
        if environment.enzyme_type != 'RNAP':
            continue

        # If we ran out of the enzyme in the environment, then there's nothing to do
        if environment.concentration < 0.0:
            continue

        # Go through sites
        for j, site in enumerate(environment.site_list):
            # We will do the binding process in this order:
            # Check site and model to use.
            # And calculate binding probability and if it'll bind
            # Check if binding site is available
            # If there are multiple enzymes that want to bind but their ranges overlap, we must choose
            # TODO: In the future, check how to decide between overlapping sites. - with a list of ranges,
            #  probabilities and enzymes

            # TODO: Generalize it for other enzymes like lacs, and topos
            # For now, only genes!
            # -----------------------------------------------------------
            if site.site_type != 'gene':
                continue

            # Get superhelical density at site
            enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= site.start][-1]
            site_superhelical = enzyme_before.superhelical

            # According model, calculate the binding probability
            # TODO: for checking the binding of overlapping enzymes, maybe I can create a dictionary of site.
            #  Or maybe with a list of ranges, sizes and probabilities, and then decide.
            #  Also figure a way how to incorporate the concentration in the enviroment and all that
            # -----------------------------------------------------------
            # We need to figure out 1 or 2 models.
            # A model for the rate in case it is modulated by supercoiling
            # And a model for calculating the binding probability.

            # Simple poisson process (constant binding)
            if site.site_model == 'poisson' or site.site_model == 'Poisson':
                binding_probability = P_binding_Poisson(site.k_min, dt)

            # TODO: Check how the topo's are going to figure this out, and how to include the environment.
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
            elif site.site_model == 'none' or site.site_model == 'None' or site.site_model == None:
                continue
            else:  # If there's no model, there's no binding
                continue

            # Check if site is available
            # -------------------------------------------------------------
            site_available = check_site_availability(site, enzyme_list, environment.size)
            if not site_available:
                continue

            # Decide if the enzyme will bind
            # -------------------------------------------------------------
            urandom = np.random.uniform()  # we need a random number

            if urandom <= binding_probability:  # and decide

                # Add enzyme
                # --------------------------------------------------------
                # We first need to figure out the position, twist and superhelical (these last two will be sorted in
                # the circuit module

                # position:
                if site.direction > 0:
                    position = site.start
                else:
                    position = site.start - environment.size

                # Create enzyme, and note that it is missing twist and the superhelical density.
                # I think it's better to fix it in the circuit module
                enzyme = Enzyme(e_type=environment.enzyme_type, name=environment.name, site=site, position=position,
                                size=environment.size, twist=0.0, superhelical=0.0)

                new_enzymes.append(enzyme)

    return new_enzymes


# ----------------------------------------------------------
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
    if site.direction > 0:
        my_range = [site.start, site.start + size]
    else:
        my_range = [site.start, site.start - size]

    # If any of them intersect
    if (set(range_before) & set(my_range)) or (set(range_after) & set(my_range)):
        available = False
    # there is an intersection
    else:
        available = True

    return available


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
    U = 1.0 / (1.0 + np.exp((sigma - SM_sigma_t) / SM_epsilon_t))  # the energy required for melting
    f = np.exp(SM_m * U)  # the activation curve
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
    U = opening_energy(sigma, *opening_p)  # the energy required for melting
    U0 = opening_energy(sigma0, *opening_p)  # energy for melting at reference sigma0
    # (should be the sigma at which k0=basal_rate
    #  was measured...)
    DU = U - U0  # Energy difference
    f = np.exp(-DU / (EE_alpha))  # the activation curve
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

