import numpy as np
import params

#---------------------------------------------------------------------------------------------------------------------
#DESCRIPTION
#---------------------------------------------------------------------------------------------------------------------
#This module contains the mathematical functions that compose the statistical part of my model
#It comprasis the necessary equations required for simulating the stochastic binding, and the
#topoisomerases/gyrases activities

#All parameters are already in the params module, but I prefer to have them here with more simple names:
v0     = params.v0
w0     = params.w0
gamma  = params.gamma
dt     = params.dt

topo_w = params.topo_w
topo_t = params.topo_t
topo_c = params.topo_c
topo_k = params.topo_k

gyra_w = params.gyra_w
gyra_t = params.gyra_t
gyra_c = params.gyra_c
gyra_k = params.gyra_k

kBT = 310.0 * params.kB_kcalmolK #The Boltzmann constant multiplied by 310K which is the temperature
                               #at which the SIDD code is ran...

#Sam Meyer's PROMOTER CURVE (parameters taken from Houdaigi NAR 2019)
SM_sigma_t   = params.SM_sigma_t
SM_epsilon_t = params.SM_epsilon_t
SM_m         = params.SM_m

#EFFECTIVE ENERGY PROMOTER CURVE (inspired by Houdaigi NAR 2019)
EE_alpha = params.EE_alpha


#---------------------------------------------------------------------------------------------------------------------
#BINDING FUNCTIONS
#---------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------
#This function is in charge of administrating the 
#probabilistic model we will use for binding events.
#model = poisson or not (non homogeneus Poisson)
#genome = the dataframe with all the information in the piece of "genome"
#sigma  = supercoiling density that the binding site
#f_ini  = initiation function
def P_binding(genome, sigma):
    f_ini = genome['fun']
    model = genome['model']
    rates = [genome['k_min'],genome['k_max']]
    opening_p = genome['params']
    if model == 'poisson':
        probability = P_binding_Poisson(rates[0])
    else:
        probability = P_binding_Nonh_Poisson( rates, sigma, f_ini, *opening_p)
    return probability

#----------------------------------------------------------
#This equation calculates the probability of binding according
#the Poisson process
def P_binding_Poisson(rate):

    rdt = rate*dt #it is what is in the exponent (is that how you call it?)
    probability = rdt*np.exp(-rdt)

    return probability

#----------------------------------------------------------
#This equation calculates the probability of binding according
#a nonhomogeneouss Poisson process, which is basically a Poisson process
#with variable rate (simply modelling).
#basal_rate - array with [0] = min_rate & [1] = max_rate
#sigma - supercoiling density
#f_ini - function of 
def P_binding_Nonh_Poisson(basal_rate, sigma, f_ini, *opening_p):

    rate = f_ini( basal_rate, sigma, *opening_p ) #the basal rate or better said, 
                                                  #the measured rate is modulated
                                                  #through a promoter activation curve (sigmoid)

    probability = rate*dt                         #The smaller dt the more accurate it is. 

    return probability

#----------------------------------------------------------
#The promoter activation curve according Sam Meyer 2019
#basal_rate[0] = minimum rate
#basal_rate[1] = maximum rate
#For this function, we use the minimum rate
def promoter_curve_Meyer( basal_rate, sigma, *opening_p ):

    U = 1.0/ ( 1.0 + np.exp( (sigma- SM_sigma_t)/SM_epsilon_t ) )  #the energy required for melting
    f = np.exp( SM_m*U )                                        #the activation curve
    rate = basal_rate[0]*f                                      #and the rate modulated through the activation curve
    return rate

#----------------------------------------------------------
#The supercoiling dependant opening energy of the promoter
#sequence. It follows a sigmoid curve, and should be
#calculated with the SIDD algorithm and following the methods
#of Houdaigui 2021 for isolating the discriminator sequence.

#Parameters:
# sigma - supercoiling density
# a,b,sigma_t,epsilon - sigmoid curve fitted parameters
def opening_energy(x, a, b, sigma_t, epsilon):
    return a+b/(1 + np.exp( -(x - sigma_t)/epsilon ) )

#----------------------------------------------------------
#The promoter activation curve relaying on the effective
#thermal energy. This curve is parametrized by the fitting
#of the openning energy, and inspired by according Sam Meyer 2019
#basal_rate[0] = minimum rate
#basal_rate[1] = maximum rate
#For this function, we use the minimum rate
def promoter_curve_opening_E( basal_rate, sigma, sigma0, *opening_p ):

    U  = opening_energy( sigma,  *opening_p ) #the energy required for melting
    U0 = opening_energy( sigma0, *opening_p ) #energy for melting at reference sigma0 
                                              # (should be the sigma at which k0=basal_rate 
                                              #  was measured...)
    DU = U-U0                        #Energy difference        
    f = np.exp( -DU/(EE_alpha) )     #the activation curve
    rate = basal_rate[0]*f
    return rate

#----------------------------------------------------------
#The promoter activation curve parametrized by the 
#opening energy fitted parameters, and the observed 
#maximum and minimum rates.
#opening_p - opening energy fitted parameters
#basal_rate[0] = minimum rate
#basal_rate[1] = maximum rate
def promoter_curve_opening_E_maxmin( basal_rate, sigma, *opening_p ):

    a = np.log( basal_rate[0]/basal_rate[1] )
   # a = np.log( max_rate/min_rate )
    b = 1+np.exp( -(sigma-opening_p[2])/opening_p[3] )
    rate = basal_rate[1] * np.exp(a/b)
    #rate = min_rate * np.exp(a/b)
    return rate

#----------------------------------------------------------
#Basically, is the same function as the previous one,
#but this one has the sigmoid inverted, hence, we
#need to adjust the location of the minimum/maximum rates in
#the equation.
#opening_p - opening energy fitted parameters
#basal_rate[0] = minimum rate
#basal_rate[1] = maximum rate
def promoter_curve_opening_E_maxmin_i( basal_rate, sigma, *opening_p ):

    a = np.log( basal_rate[1]/basal_rate[0] )
    b = 1+np.exp( -(sigma-opening_p[2])/opening_p[3] )
    rate = basal_rate[0] * np.exp(a/b)
    return rate


#---------------------------------------------------------------------------------------------------------------------
#TOPOISOMERASE ACTIVITY FUNCTIONS
#---------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------
#Calculates the amount of coils removed by topoisomerase I
#activity. This function only depends on the supercoiling density (sigma)
#I took this function from Sam Meyer's paper (2019)
def topoisomerase_activity(sigma):
    #the function has the form of (concentration*sigmoid)*rate*dt
    a = topo_c*topo_k*dt
    b = 1+np.exp( (sigma - topo_t)/topo_w )
    sigma_removed = a/b
    return sigma_removed

#----------------------------------------------------------
#Calculates the amount of coils removed by gyrase
#activity. This function only depends on the supercoiling density (sigma)
#I took this function from Sam Meyer's paper (2019)
def gyrase_activity(sigma):
    #the function has the form of (concentration*sigmoid)*rate*dt
    a = gyra_c*gyra_k*dt
    b = 1+np.exp( -(sigma - gyra_t)/gyra_w )
    sigma_removed = -a/b
    return sigma_removed

