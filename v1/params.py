import numpy as np

#---------------------------------------------------------------------------------------------------------------------
#DESCRIPTION
#---------------------------------------------------------------------------------------------------------------------
#This module contains all the parameters that will be used in the main code
#---------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
#PARAMETERS
#---------------------------------------------------------------------------------------------------------------------
#Physical constants
#kB_kcalmolK = 1/310
kB_kcalmolK = 1.987204259 * .001 #10**(-3) #Boltzman constant in kcal/(mol*K) units...
#kB_kcalmolK = .593 #10**(-3) #Boltzman constant in kcal/(mol*K) units...
dt     = 1.0            #timestep (seconds)
v0     = 30.0           #Velocity (bp/sec) of RNAPs
T_bp   = 10.5           #Number of bp per helical turn
w0     = 2.0*np.pi/T_bp #Relaxed twist density per bp
gamma  = 0.2*w0         #How much supercoiling is inyected per bp
#sigma0 = -0.06          #Initial supercoiling density

#TOPOISOMERASE I
topo_w = 0.012#width
topo_t = -0.04#thresholds
topo_c = 0.25#0.1#concentration micromolar 0.025 in meyer -> this is too negative...
topo_k = 0.001 #basal rate

#GYRASE
gyra_w = 0.025 #width
gyra_t = 0.01 #threshold
gyra_c = 0.25#.01 #concentration micromolarb 0.25 in meyer
gyra_k = 0.001 #minus because it removes negative supercoils
              #I think it makes more sense to put the negative in the equation rather than in
              #the parameter

#Sam Meyer's PROMOTER CURVE (parameters taken from Houdaigi NAR 2019)
SM_sigma_t   =  -0.042 #threshold of promoter openning
SM_epsilon_t =   0.005 #width of the crossover
SM_m         =   2.5   #effective thermal energy that sets the SC activation factor

#EFFECTIVE ENERGY PROMOTER CURVE (inspired by Houdaigi NAR 2019)
EE_alpha = 3.3 #The efective energy is beta = kBT/alpha...
               #In the equation we use kBT is canceled out, that's why I only use the factor... 



#---------------------------------------------------------------------------------------------------------------------
#OBJECTS (proteins) paramters
#---------------------------------------------------------------------------------------------------------------------

#OBJECT SIZES - You cand add/modify NAPs or RNAPolymerase. 
OBJECT_size  = {"ori":30, "lacI":20, "IHF":30, "FIS":10, "RNAP":30, "EXT_L":0, "EXT_R":0}
#I assume lacI is as big as its binding site.
#I assume ori and RNAP are as big as promoter sequences ~ 30bp
#EXT - It is a fake protein that I add to simulate the boundaries
#ori - origin of replication is treated as a NAP as it is treated as barrier because proteins bind there and 
#      it is possible that it is anchored to the membrane

