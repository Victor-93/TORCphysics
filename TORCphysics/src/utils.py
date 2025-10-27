import numpy as np
import pandas as pd
from TORCphysics import params
import sys

# TODO: Decide which of these parameters you need
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

def calculate_twist_writhe_partition(superhelical, length, twist_writhe_ratio=params.twist_writhe_ratio):
    """
    Calculate the twist–writhe partition for a given superhelical density.

    This function computes the linking number difference (ΔLk) and partitions it
    into twist and writhe contributions according to a given twist–writhe ratio.
    The calculation can be applied to either local or global superhelical densities.

    Parameters
    ----------
    superhelical : float
        Local or global superhelical density (σ).
    length : float
        Length of the topological domain in base pairs (bp).
    twist_writhe_ratio : float, optional
        Fraction of the linking difference assigned to twist (default from `params.twist_writhe_ratio`).

    Returns
    -------
    dLk : float
        Linking number difference (ΔLk).
    twist : float
        Partitioned twist contribution (in turns).
    writhe : float
        Partitioned writhe contribution (in turns).
    """
    L0 = length / params.T_bp  # Relaxed linking number for a B-DNA structure of same length
    dLk = superhelical * L0    # Linking difference
    twist = dLk * twist_writhe_ratio
    writhe = dLk - twist
    return dLk, twist, writhe

# Having an enzyme E sourrounded by a molecules L and R, on the left and right respectively (L____E____R),
# calculate the change in twist on the left and right sides of E given the fact that enzyme E completely
# leaks the twist on both sides (it doesn't form a topological barrier).
def instant_twist_transfer(z,z_list):

    # Get enzyme before z (z_b, on the left) and enzyme after z (z_a, on the right)
    z_b = get_enzyme_before_position(position=z.position - 10, enzyme_list=z_list)
    z_a = get_enzyme_after_position(position=z.position + 10, enzyme_list=z_list)

    # Total twist trapped in the region
    total_twist = z_b.twist + z.twist
    # Calculate lengths
    length_left = calculate_length(z_b, z)
    length_right = calculate_length(z, z_a)

    # Partitionate the twist
    #  *** Small update. It is very unlikely but sometimes it can divide by zero. Let's ignore the cases when
    #  this happens and don't do anything
    b = length_left + length_right
    if b == 0:
        return 0, 0

    twist_left = total_twist * length_left / b #(length_left + length_right)
    twist_right = total_twist * length_right / b# (length_left + length_right)

    # And calculate the actual change in twist
    dtwist_left = twist_left - z_b.twist
    dtwist_right = twist_right - z.twist

    return dtwist_left, dtwist_right


# Having RNAP molecule z and a list of enzymes z_list on the DNA, it calculates the updated position and the
# twist induced on both sides of z (left and right) according a mechanical model where the RNAP can stall
# according the balance of torques.
def RNAP_stall_mec_model(z, z_list, dt):
    # Get enzymes on the right and left
    z_right = [e for e in z_list if e.position > z.position][0]  # after - On the right
    z_left = [e for e in z_list if e.position < z.position][-1]  # before - On the left

    # First, we need to calculate the Torque acting on our RNAP.
    # Calculate torques and determine if the RNAP will stall
    torque_right = Marko_torque(z.superhelical)  # Torque on the right
    torque_left = Marko_torque(z_left.superhelical)  # Torque on the left

    if z.direction > 0:  # Moving to the right
        torque = torque_right - torque_left
    else:  # Moving to the left
        torque = torque_left - torque_right

    velocity = velocity_2022SevierBioJ(z, torque)

    position = z.direction * velocity * dt

    # Injects twist: denatures w = gamma*velocity*dt base-pairs
    twist_left = -z.direction * z.effect_model.gamma * params.w0 * velocity * dt
    twist_right = z.direction * z.effect_model.gamma * params.w0 * velocity * dt

    # Check if there's one enzyme on the direction of movement. If there is one, then it will stall to avoid
    # clashing
    if z.direction > 0:  # Moving to the right
        if z_right.position - (z.position + position) <= 0:
            position = 0.0
            twist_left = 0.0
            twist_right = 0.0
    else:  # Moves to the left
        if z_left.position - (z.position + position) >= 0:
            position = 0.0
            twist_left = 0.0
            twist_right = 0.0

    return position, twist_left, twist_right

# Having the case of E1___E____E2, being the E's the enzymes in a region, calculate the twist and superhelical
# density surrounding enzyme E. This function is useful for calculating the effects of topoisomerase or enzymes
# that act on regions.
def get_superhelical_in_region(z, z_list):
    z_b = get_enzyme_before_position(position=z.position - 10, enzyme_list=z_list)  # Get the enzyme before
    z_a = get_enzyme_after_position(position=z.position + 10, enzyme_list=z_list)  # Get the enzyme after z
    total_twist = z_b.twist + z.twist  # Total twist in the region.
    # This is the total superhelical density of a region; enzyme X does not block supercoils  O______X_____O
    superhelical = total_twist / (params.w0 * (z_a.position - z_b.position))
    return superhelical, total_twist


# Sidmoid function used for modulating the rate in the thermodynamic model of strand separation.
# This function is taken from Houdagui et al. 2019
# The parameter SM_m is an effective thermal energy (1/SM_m). Basically, it used to limit the negative exponential
# to values from 0 to 1, so the rate can follow the same shape as the free energy of melting.
def opening_function(superhelical, threshold, width):
    return params.SM_m / (1.0 + np.exp(-(superhelical - threshold) / width))

# TODO: Test and document
# n = spacer length ( dimension less)
# superhelical = local superhelical density
# twist_o = optimal twist angle for RNAP binding (degrees)
# NOTES: 1.- averaged B-DNA twist angle (twist_deg) is also in degrees.
#        2.- Output is in kBT units
# TODO: This equation and parameters are taken from the NAR paper Forquet, R., et al. (2022).
#  But the units are a bit confusin. the twist angle are in degrees, but the stiffness constant in radians.
def spacer_length_free_energy(n, superhelical, twist_o):
    a = twist_o/n - params.twist_deg * (1.0+superhelical)
    spacer_energy = 0.5 * n * params.k_twist_kBT * a**2
    return spacer_energy


# TODO: Test and document
# Calculates the energy spent of going from state 1 with superhelical 1 to state 2 with superhelical2.
def change_supercoiling_free_energy(superhelical1, superhelical2, length):
    a = superhelical2 * superhelical2 - superhelical1 * superhelical1
    return params.q * a / (2 * length)


# TODO: Test and document
# Calculates superhelical free energy
def superhelical_free_energy(superhelical, length):
    return params.q * superhelical * superhelical / (2 * length)


# The idea of this module, is to define utility functions, that are used several times by the other main modules.
# TODO: You need to test this!
def get_enzyme_before_position(position, enzyme_list):
    #    enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= position][-1]
    enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= position][-1]

    # I did this to test, and sometimes for circular structures this can produce errors. But maybe is not due to
    # this function, it might be more of the circular nature. Anyway, we also need to test this function independently
    # Of write some warnings.
    #n = len([enzyme for enzyme in enzyme_list if enzyme.position <= position])
    #if n >0:
    #    enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= position][-1]
    #else:
    #    enzyme_before = None
    #    a=2
    return enzyme_before


def get_enzyme_after_position(position, enzyme_list):
    enzyme_after = [enzyme for enzyme in enzyme_list if enzyme.position >= position][0]
    return enzyme_after


# site = site to be bound
# environmental = enzyme binding the site
# enzyme_list = list of enzymes currently bound to the DNA
def check_site_availability(site, environmental, enzyme_list):
    # Enzyme before and after the start site
    enzyme_before = get_enzyme_before_position(position=site.start, enzyme_list=enzyme_list)
    enzyme_after = get_enzyme_after_position(position=site.start, enzyme_list=enzyme_list)

    # Let's build ranges
    # ----------------------------------------------------
    # Range of the enzyme to bind
    bind_a, bind_b = get_enzyme_to_bind_ranges(site=site, environmental=environmental)
    # Range of the enzyme on the left (before the start site)
    before_a, before_b = get_enzyme_ranges(enzyme=enzyme_before)
    # Range of the enzyme on the right (after the start site)
    after_a, after_b = get_enzyme_ranges(enzyme=enzyme_after)

    # Check if ranges overlap
    # ----------------------------------------------------
    # Let's check if it overlaps with the enzyme before the start site (left)
    overlap = max(before_a, bind_a) <= min(before_b, bind_b)
    if overlap:
        available = False
        return available

    # Let's check if it overlaps with the enzyme after the start site (right)
    overlap = max(bind_a, after_a) <= min(bind_b, after_b)
    if overlap:
        available = False
        return available

    # If the code reach until here, then it means that the ranges didn't overlap, so the site is available!
    available = True
    return available


# This calculates the ranges of enzymes. Each enzyme has the range (a,b). This range represents the amount of space
# that the enzymes occupy the DNA. It is not the space that the enzymes actually occupy...
def get_enzyme_ranges(enzyme):
    a = enzyme.position - (enzyme.size - enzyme.effective_size) * 0.5
    b = enzyme.position + (enzyme.size + enzyme.effective_size) * 0.5
    return a, b


# Calculates the range that will cover an environmental trying to bind a site. This range is in the form (a,b).
def get_enzyme_to_bind_ranges(site, environmental):
    if site.direction == 1:
        a = site.start - (environmental.size + environmental.effective_size) * 0.5
        b = site.start + (environmental.size - environmental.effective_size) * 0.5
    elif site.direction == 0 or site.direction == -1:
        a = site.start - (environmental.size - environmental.effective_size) * 0.5
        b = site.start + (environmental.size + environmental.effective_size) * 0.5
    else:
        raise ValueError('Error, invalid direction in site ' + site.name)
    return a, b


# Direction ==1
# Enzymes that advance to the right (RNAPs with ->>> direction of transcription) load/bind just behind the start site:
#  ___|RNAP|START|_______ ==== ____|EffectiveSize|START|________, Notice that the actual size (not effective size) might
#  overlap with the START size, but the actual contact happens just before the start site
# Direction == 0,-1
# Enzymes that don't move (do not have direction) or that advance to the left (RNAP with <<<- direction of
# transcription) load/bind just after the start site:
# ___|START|RNAP________ === ___|START|EffectiveSize|_____; Notice that the complete size could overlap with the start
# site, but the contact with the DNA happens just at the start of the site.
def get_enzyme_to_bind_position(site, environmental):
    if site.direction == 1:
        position = site.start - environmental.effective_size
    elif site.direction == 0 or site.direction == -1:
        position = site.start
    else:
        raise ValueError('Error, invalid direction in site ' + site.name)
    return position


def new_enzyme_start_position(site, environmental):
    if site.direction == 1:
        position = site.start - environmental.effective_size
    elif site.direction == 0 or site.direction == -1:
        position = site.start
    else:
        raise ValueError('Error, invalid direction in site ' + site.name)
    return position


# ----------------------------------------------------------
# This function calculates the length between two objects (proteins) considering their effective size.
# Basically, according the effective size is the size of the enzyme that actually touches the DNA.
def calculate_length(z0, z1):
    x0 = z0.position  # positions
    x1 = z1.position
    b0 = z0.effective_size
    #    b0 = z0.size  # size -_-
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
# to the length between object Z0 and Z1.
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
    # b_n = zn.size
    b_n = zn.effective_size
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


# This equation calculates the probability of binding according
# a Non-homogeneous Poisson process, which is basically a Poisson process
# with variable rate (simply modelling).
# It assumes that the rate was already varied and obtained by one of the opening energies
# sigma - supercoiling density
def P_binding_Nonh_Poisson(rate, dt):
    probability = rate * dt  # The smaller dt the more accurate it is.

    return probability


def Poisson_process(rate, dt):
    """
    Calculates probability of a Poisson process. Note that this is an approximation for a timestep dt that is smaller
    than the rate. Hence, it calculates the probability of observing one occurrence.

    Parameters
    ----------
    rate : float
        This is the frequency (rate) at which one event occurs (1/s).
    dt : float
        Timestep in seconds (s).

    Returns
    ----------
    probability : float
        It represents the probability of observing one occurrence.
    """
    rdt = rate * dt  # it is what is in the exponent (is that how you call it?)
    probability = rdt * np.exp(-rdt)
    return probability


# ----------------------------------------------------------
# The promoter activation curve according Sam Meyer 2019
# For this function, we use the minimum rate
def promoter_curve_Meyer(basal_rate, superhelical):
    u = 1.0 / (1.0 + np.exp((superhelical - SM_sigma_t) / SM_epsilon_t))  # the energy required for melting
    f = np.exp(SM_m * u)  # the activation curve
    rate = basal_rate * f  # and the rate modulated through the activation curve
    return rate


def read_csv_to_dict(filename):
    """
    Reads csv file and puts it in a dictionary
    """
    return pd.read_csv(filename).to_dict()


def site_match_by_name(site_list, label):
    """ Given the site_list, filters sites by name 'label'.

    Parameters
    ----------
    site_list : list
        It is a list of Sites.
    label : str
        Name of site the enzyme is bound to.

    Returns
    ----------
    list : The site with the name 'label'.

    """

    if label in [site.name for site in site_list]:
        for site in site_list:
            if site.name == label:
                return site  # the first one?
    else:
        return None


def site_match_by_type(site_list, label):
    """ Given the site_list, filters sites by site_type 'label'.

    Parameters
    ----------
    site_list : list
        It is a list of Sites.
    label : str
        Type of site.

    Returns
    ----------
    list : A list of sites of the type 'label'.

    """
    #        enzyme_before = [enzyme.position for enzyme in enzyme_list if enzyme.position <= site.start][-1]
    site_list = [site for site in site_list if site.site_type == label]
    return site_list


# Read fasta file. Returns the sequence
def read_fasta(file_name):
    fasta_file = open(file_name, 'r')
    header = fasta_file.readline()  # Reads the header
    lines = []  # This one will contain all the lines
    while True:
        line = fasta_file.readline()
        if not line:  # If we reach the end, break loop
            break
        lines.append(line[:-1])  # -1 so we remove the spacing
    sequence = ''.join(lines)  # And join all lines to obtain sequence
    fasta_file.close()
    return sequence

# Torque calculated using Marko's elasticity model
def Marko_torque(sigma):
    if abs(sigma) <= abs(params.sigma_s):
        torque = sigma * params.cs_energy / params.w0
    elif abs(params.sigma_s) < abs(sigma) < abs(params.sigma_p):
        torque = np.sign(sigma)*np.sqrt(
            2 * params.p_stiffness * params.g_energy / (1 - params.p_stiffness / params.cs_energy)) / params.w0_nm
        # torque = np.sqrt(
        #     2 * params.p_stiffness * params.g_energy / (1 - params.p_stiffness / params.cs_energy)) / params.w0_nm
    elif abs(sigma) >= abs(params.sigma_p):
        torque = sigma * params.p_stiffness / params.w0
    else:
        #print('Error in Marko_torque function')
        torque = 0.0
        # sys.exit()
    return torque


#  The velocity has the form: v = vmax/ (1+e^{k(T_0 - T_c)} )
#  where vmax = maximum velocity, k = torque parameter, T_0 = Torque acting on enzyme
#  and T_c = cutoff or stalling torque.
#  This function is based on the 2022SevierBioJ paper
def velocity_2022SevierBioJ(z, input_torque):
    # top = 2.0 * z.effect_model.velocity
    #torque = input_torque
    torque = abs(input_torque) # This is to make it symmetric. We do it like this to avoid cases with huge torques
    top = z.effect_model.velocity
    exp_arg = z.effect_model.kappa * (torque - z.effect_model.stall_torque)
    exp_arg = np.float128(exp_arg)

    # Define a maximum value for the argument to exp to prevent overflow
    max_exp_arg = 709  # slightly below the overflow threshold for float64
    # Clip the argument to the maximum value
    exp_arg_clipped = np.clip(exp_arg, a_min=None, a_max=max_exp_arg)
    #down = 1.0 + np.exp(exp_arg)
    down = 1.0 + np.exp(exp_arg_clipped)
    velocity = top / down
    return velocity

