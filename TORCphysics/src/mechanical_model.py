import numpy as np
import pandas as pd
import params
import sys

#---------------------------------------------------------------------------------------------------------------------
#DESCRIPTION
#---------------------------------------------------------------------------------------------------------------------
#This module contains the mathematical functions that compose the mechanical part of my model
# It also includes the functions that update the fake ends, as well as how the addition of
# RNAPs modify the topology of our system.

#All parameters are already in the params module, but I prefer to have them here with more simple names:
v0     = params.v0
w0     = params.w0
gamma  = params.gamma
dt     = params.dt


#----------------------------------------------------------
#This function calculates the length between two objects (proteins) considering their sizes
def calculate_length(Z0, Z1):
    x0 = Z0['pos']  #positions
    x1 = Z1['pos']
    b0 = Z0['size'] #size -_-
    b1 = Z1['size']
    #There are 4 posibilities
    if Z0['direction'] >= 0 and Z1['direction'] >= 0:
        length = (x1 - b1) - x0
    elif Z0['direction'] >= 0 and Z1['direction'] < 0:
        length = x1 - x0
    elif Z0['direction'] < 0 and Z1['direction'] >= 0:
        length = (x1 - b1) - (x0 + b0)
    elif Z0['direction'] < 0 and Z1['direction'] < 0:
        length = x1 - (x0 + b0)
    else:
        print("Something went wrong in lengths")
        sys.exit()
    #length+=1 #1 bp needs to be added
    return length

#----------------------------------------------------------
#This function calculates/updates the twist parameter according
#the supercoiling value of the current object Z0, and according
#the length between object Z0 and Z1.
def calculate_twist(Z0,Z1):
    length = calculate_length(Z0,Z1) #First, I need to get the length
    sigma = Z0['superhelical']
    twist  = sigma*w0*length
    return twist

#----------------------------------------------------------
#This function calculates/updates the supercoiling according
#the twist of the current object Z0, and the distance between
#Z1-Z0
def calculate_supercoiling(Z0, Z1):
    length = calculate_length(Z0,Z1) #First, I need to get the length
    twist  = Z0['twist']             #and twist
    if length != 0:
        sigma  = twist/(w0*length)       #and calculate the supercoiling
    else:
        sigma = 0 #I don't know if this is a solution... #But basically, this happens when a RNAP
                  #binds right after one has bound
    #if length == 0:
    #    print("Aqui")
    #    print(Z0)
    #    print(Z1)
    #    sys.exit()
    return sigma

#----------------------------------------------------------
#This function calculates/updates the position of an RNAP Z
def RNAP_motion(Z):

    Zn = Z['start'] + Z['direction']*v0*dt  #Zn = new Z (new position)

    return Zn

#----------------------------------------------------------
#This function calculates twist inyected on the left and right
#by the RNAP Z
def twist_injected(Z):

    twist_left  = -Z['direction']*gamma*v0*dt
    twist_right =  Z['direction']*gamma*v0*dt

    return twist_left, twist_right


#-----------------------------------------------------------------------
#Get's the start and end positions of the fake boundaries (for circular DNA) 
# In case that there is not fake boundaries, Z_N should be the last element [-1],
# in case that you have N objects including the fake boundaries, Z_N -> [N-2]
def get_start_end_c(Z_0, Z_N, nbp):

    b_0 = Z_0['size']
    b_N = Z_N['size']
    x_0 = Z_0['pos']                        #position of first object
    x_N = Z_N['pos']                        #position of last object

    #fake start
    if  Z_N['direction'] >= 0: #depends on the direction
        start_0 = 0 - (nbp - x_N)       #this is the start site of the fake bit,
    else:
        start_0   = 0 - ( nbp - (x_N + b_N) ) # the size of the last object is considered

    #fake end
    if  Z_0['direction'] >= 0: #depends on the direction
        start_N   = nbp + x_0 - b_0     #I think I had the sign wrong...
    else:
        start_N   = nbp + x_0

    return start_0, start_N

#-----------------------------------------------------------------------



