import numpy as np
import pandas as pd
from TORCphysics import params

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


# The idea of this module, is to define utility functions, that are used several times by the other main modules.
# TODO: You need to test this!
def get_enzyme_before_position(position, enzyme_list):
    #    enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= position][-1]
    enzyme_before = [enzyme for enzyme in enzyme_list if enzyme.position <= position][-1]
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
