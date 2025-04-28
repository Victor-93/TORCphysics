import numpy as np

# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This module contains all the parameters that will be used in the main code

# ---------------------------------------------------------------------------------------------------------------------
# TODO: ORGANIZE PARAMETERS by
#  Physical Constants
#  Binding Model params
#  Effect Models params
#  Unbinding Models params
#  Plotting constants
# ---------------------------------------------------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------

# Physical constants
# -----------------------------------------------------------
Temperature = 300  # K - temperature
# kB_kcalmolK = 1/310
kBT_pN_nm = 4.1  # pN nm at T=300K
kB_kcalmolK = 1.987204259 * .001  # 10**(-3) #Boltzman constant in kcal/(mol*K) units...
# kB_kcalmolK = .593 #10**(-3) #Boltzman constant in kcal/(mol*K) units...
dt = 1.0  # Default timestep (seconds)
bp_nm = .34  # nm - base-pair rise
T_bp = 10.5  # Number of bp per helical turn
w0 = 2.0 * np.pi / T_bp  # Relaxed twist density per bp (rad/bp)
w0_nm = 2.0 * np.pi / (T_bp * bp_nm)  # rad * nm^{-1} - Now in distance units
R = 1.9872  # Gas constant cal⋅K−1⋅mol−1
RT = R*Temperature/1000.0  # kcal * mol^{-1}
q = 2350 * RT  # Coefficient associated to the superhelical free energy
#ATP_hydrolysis = 28  # kJ / mol
ATP_hydrolysis = 6.7  # kcal / mol

# Elasticity parameters - from Marko's elasticity model
# -----------------------------------------------------------
P_length = 24.0  # nm - twist stiffness of plectonomic DNA - is it a persistence length?
C_length = 95.0#100.0  # nm - twist persistence length - 95nm used in Marko's paper
A_length = 50.0  # nm - persistence length
f_stretching = 1.0  # 0.15  # pN - stretching force of enzymes in vivo - this might not be correct, but it is in the range of
# low forces in the Marko's elasticity model
p_stiffness = kBT_pN_nm * P_length * w0_nm * w0_nm  # pN - stiffness related to P
c_stiffness = kBT_pN_nm * C_length * w0_nm * w0_nm  # pN - stiffness related to C
# related to free energy of twist stiffness of extended state
cs_energy = c_stiffness * (
        1 - ((C_length / 4 * A_length) *
             np.sqrt(kBT_pN_nm / (A_length * f_stretching))))
# related to free energy of stretched state
g_energy = f_stretching - np.sqrt(kBT_pN_nm * f_stretching / A_length)
# |sigma| <= |sigma_s| - > only twist exists
sigma_s = (1 / cs_energy) * np.sqrt(2 * p_stiffness * g_energy / (1 - p_stiffness / cs_energy))
# |sigma| > |sigma_p| - > only writhe exists
sigma_p = (1 / p_stiffness) * np.sqrt(2 * p_stiffness * g_energy / (1 - p_stiffness / cs_energy))

# From Forquet, R., et al. (2022). Quantitative contribution of the spacer length in the supercoiling-sensitivity of bacterial promoters. Nucleic Acids Research, 50(13), 7287–7297. https://doi.org/10.1093/nar/gkac579
k_twist_kBT = 71.4 # kBT *rad^{-2} - In the paper it says that it is rads, but shouldn't it be degrees? I'll have to check
                   # This expression is in kBT units
k_twist_pN_nm = 71.4 * kBT_pN_nm # And now multiplied by kBT but in pN nm units.
twist_deg = 34.3 # B-DNA twist in degrees
spacer_factor = 15.0 # Effective energy used to scale down the geometric modulation of binding rate according the
                     # spacer length. We chose this factor as it limits the binding at hyper superhelical values.

# ---------------------------------------------------------------------------------------------------------------------
# SPECIFIC ENZYME PARAMETERS
# ---------------------------------------------------------------------------------------------------------------------
# RNA Polymerase (RNAP)
v0 = 30.0  # 60.0  # Velocity (bp/sec) of RNAPs
gamma = 0.1#0.157#0.13#0.05 #0.835 #0.5 # How much supercoiling is injected per bp
stall_torque = 12.0#12.0  # 10.5 * 5  # * 17 # pN * nm - from Gleng's papers which cited another paper.
# 12pN*nm according 2022SevierBioJ
sigma_stall = 0.6  # If sigma greater than this, then the RNAP will stall - According Gleng?
RNAP_kappa = 0.5  # 12pN^{-1} - According 2022SevierBioJ. This is a parameter used in calculating torque dependant
# velocity.

# TOPOISOMERASE I
topo_b_k_off = 0.5

# GYRASE
gyra_b_k_off = 0.5

# Thresholds
# ---------------------------------------------------------------------------------------------------------------------
gene_RNAP_threshold = (60/12) * v0 # This quantity is used to define a threshold in the number of bound RNAPs to a gene.
                                   # It takes the form: (60secs/ ntranscripts) * v0, which gives a distance that
                                   # indicates the maximum number of RNAPs per gene section that can be bound.
                                   # Example: 60/12 * 30 = 150 bp, that means that the threshold for number of bound
                                   # RNAPs is gene_length/150; which for a gene_length=900 is 6.

# ---------------------------------------------------------------------------------------------------------------------
# STOCHASTIC TOPOISOMERASE  PARAMETERS
# ---------------------------------------------------------------------------------------------------------------------
# TopoIRecognition
topo_b_w = 0.00425#0.012  # binding width
topo_b_t = 0.012##-0.04  # binding threshold
topo_b_k_on = 0.000354 #0.005

# Linear Effect
topo_e_k_cat = 15.11

# GyraseRecognition
gyra_b_w = 0.0123 #0.025  # binding width
gyra_b_t = -0.059  #0.01  # binding threshold
gyra_b_k_on = 0.000217  #0.005

# Linear Effect
gyra_e_k_cat = 13.94
gyra_e_sigma0 = -0.146

# Uniform effect
gyra_uniform_k_cat = -10.5  # (bp/second)
topoI_uniform_k_cat = 7.5  # (bp/second)

# ---------------------------------------------------------------------------------------------------------------------
# OTHER MODELS
#  ---------------------------------------------------------------------------------------------------------------------
# PoissonBinding
k_on = 0.01  # default binding rate

# PoissonUnBinding
k_off = 0.02  # default unbinding rate for any enzyme - For both topoI and gyrase, it is approximately ~0.2

# Houdagui et al. 2019 parameters for topo I activity (effect)
topo_sam_width = 0.012  # effect width
topo_sam_threshold = -0.04  # effect threshold
topo_sam_kcat = 0.001  # basal rate # k_cat
# topo_c = 0.25#0.1#concentration micromolar 0.025 in meyer -> this is too negative...


# Houdagui et al. 2019 parameters for gyrase activity (effect)
gyra_sam_width = 0.025  # effect width
gyra_sam_threshold = 0.01  # effect threshold
# gyra_c = 0.25#.01 #concentration micromolarb 0.25 in meyer - I'll use concentration in nM better
gyra_sam_kcat = 0.001  # minus because it removes negative supercoils
# I think it makes more sense to put the negative in the equation rather than in
# the parameter


# Sam Meyer's PROMOTER CURVE (parameters taken from Houdaigi NAR 2019)
SM_sigma_t = -0.042  # threshold of promoter openning
SM_epsilon_t = 0.005  # width of the crossover
SM_m = 2.5  # effective thermal energy that sets the SC activation factor

# EFFECTIVE ENERGY PROMOTER CURVE (inspired by Houdaigi NAR 2019)
EE_alpha = 3.3  # The efective energy is beta = kBT/alpha...
# In the equation we use kBT is canceled out, that's why I only use the factor...


# ---------------------------------------------------------------------------------------------------------------------
# OBJECTS (proteins) paramters
# ---------------------------------------------------------------------------------------------------------------------

# OBJECT SIZES - You cand add/modify NAPs or RNAPolymerase.
OBJECT_size = {"ori": 30, "lacI": 20, "IHF": 30, "FIS": 10, "RNAP": 30, "EXT_L": 0, "EXT_R": 0}
# I assume lacI is as big as its binding site.
# I assume ori and RNAP are as big as promoter sequences ~ 30bp
# EXT - It is a fake protein that I add to simulate the boundaries
# ori - origin of replication is treated as a NAP as it is treated as barrier because proteins bind there and
#      it is possible that it is anchored to the membrane
