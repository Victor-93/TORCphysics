import numpy as np

# ---------------------------------------------------------------------------------------------------------------------
# DESCRIPTION
# ---------------------------------------------------------------------------------------------------------------------
# This module contains all the parameters that will be used in the main code
# ---------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------------------------------------------------
# Physical constants
# -----------------------------------------------------------
# kB_kcalmolK = 1/310
kBT_pN_nm = 4.1  # pN nm at T=300K
kB_kcalmolK = 1.987204259 * .001  # 10**(-3) #Boltzman constant in kcal/(mol*K) units...
# kB_kcalmolK = .593 #10**(-3) #Boltzman constant in kcal/(mol*K) units...
dt = 1.0  # timestep (seconds)
v0 = 30.0#60.0  # Velocity (bp/sec) of RNAPs
bp_nm = .34  # nm - base-pair rise
T_bp = 10.5  # Number of bp per helical turn
w0 = 2.0 * np.pi / T_bp  # Relaxed twist density per bp
w0_nm = 2.0 * np.pi / (T_bp * bp_nm)  # nm^{-1} - Now in distance units
gamma = 0.2 * w0  # How much supercoiling is inyected per bp
# sigma0 = -0.06          #Initial supercoiling density
stall_torque = 10.5 * 5 #* 17 # pN * nm - from Gleng's papers which cited another paper
sigma_stall = 0.6  # If sigma greater than this, then the RNAP will stall

# Elasticity parameters - from Marko's elasticity model
# -----------------------------------------------------------
P_length = 24.0  # nm - twist stiffness of plectonomic DNA - is it a persistence length?
C_length = 100.0  # nm - twist persistence length
A_length = 50.0  # nm - persistence length
f_stretching = 0.15  # pN - stretching force of enzymes in vivo - this might not be correct, but it is in the range of
# low forces in the Marko's elasticity model
Temperature = 300  # K - temperature
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

k_on = 0.01  # default binding rate for any enzyme

# TOPOISOMERASE I
topo_b_w = 0.012  # binding width
topo_b_t = -0.04  # binding threshold
topo_b_k_on = 0.005
topo_b_k_off = 0.5
topo_e_w = 0.012  # effect width
topo_e_t = -0.04  # effect threshold
# topo_c = 0.25#0.1#concentration micromolar 0.025 in meyer -> this is too negative...
topo_k = 0.001 #basal rate # k_cat

# GYRASE
gyra_b_w = 0.025  # binding width
gyra_b_t = 0.01  # binding threshold
gyra_b_k_on = 0.005
gyra_b_k_off = 0.5
gyra_e_w = 0.025  # effect width
gyra_e_t = 0.01  # effect threshold
# gyra_c = 0.25#.01 #concentration micromolarb 0.25 in meyer - I'll use concentration in nM better
gyra_k = 0.001 #minus because it removes negative supercoils
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
