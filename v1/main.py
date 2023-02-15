import params
import mechanical_model as mm
import statistical_model as sm
import argparse
import numpy as np
import pandas as pd
import sys

#---------------------------------------------------------------------------------------------------------------------
#DESCRIPTION
#---------------------------------------------------------------------------------------------------------------------
#This program simulates a genetic circuit under certain conditions given the inputs:
#circuit.csv, genome.csv, objects.csv and enviroment.csv
#According these inputs, RNAPs will stochastically bind the DNA and will generate
#supercoiling accordingly.
#---------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
#PARAMETERS
#---------------------------------------------------------------------------------------------------------------------

#All parameters are already in the params module, but I prefer to have them here with more simple names:
v0     = params.v0
w0     = params.w0
gamma  = params.gamma
dt     = params.dt

posfile='positions.txt'
sigmafile='supercoiling.txt'
namefile='object.txt'


#---------------------------------------------------------------------------------------------------------------------
#INPUTS
#---------------------------------------------------------------------------------------------------------------------

# Create the parser
#---------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Version1 of the physical model oftranscription-supercoiling")
parser.add_argument("-f", "--frames", type=int, action="store", help="Number of frames (timesteps)", default=5000)
parser.add_argument("-c", "--continuation", action="store_true", help="Continuation of a simulation")
parser.add_argument("-ic", "--input_circuit", action="store", help="Circuit input file", default="circuit.csv")
parser.add_argument("-ig", "--input_genome", action="store", help="Genome input file", default="genome.csv")
parser.add_argument("-io", "--input_objects", action="store", help="Objects input file", default="objects.csv")
parser.add_argument("-ie", "--input_enviroment", action="store", help="Enviroment input file", default="enviroment.csv")
parser.add_argument("-s", "--series", action="store_true", help="Print dynamic results per timestep")
parser.add_argument("-t", "--test", action="store_true", help="Run series of tests stored in the test folder")
parser.add_argument("-o", "--output", action="store", help="Output prefix for output files", default="output")


#BORRAESTO
parser.add_argument("-l", "--linear", action="store_true", help="Include this option if linear DNA")
parser.add_argument("-im", "--initiation_model", action="store", help="Initiation method", default="maxmin")

# Process terminal commands
#---------------------------------------------------------------------------------------------------------------------
args = parser.parse_args()

#Let's put this information into variables
frames        = args.frames                   #Number of frames
c_circuit     = args.input_circuit            #Circuit file
c_genome      = args.input_genome             #Genome file
c_objects     = args.input_objects            #Objects file
c_enviroment  = args.input_enviroment         #Enviroment file
coutput       = args.output                   #Output file prefix

if args.continuation:        #If this is the continuation of a previous run
    continuation=True
else:
    continuation=False
if args.series:              #If series, then it will print dynamic output files (.txt)
    series=True
else:
    series=False
if args.test:                #If true, then we run some tests and stop
    test=True
else:
    test=False

# Read the input files (*csv) and prepare data frames
#---------------------------------------------------------------------------------------------------------------------

# Circuit
circuit_df = pd.read_csv(c_circuit,sep='\t')
# Genome
genome_df = pd.read_csv(c_genome,sep='\t')
genome_df = genome_df.drop_duplicates(subset=['name']) #Be careful with this
#genome_df = genome_df.sort_values(by=["start"])
genome_df = genome_df.reset_index(drop=True)
# Objects
o_df = pd.read_csv(c_objects,sep='\t')
o_df = o_df.sort_values(by=["pos"])
o_df = o_df.reset_index(drop=True)
# Enviroment
enviroment_df = pd.read_csv(c_enviroment,sep='\t')
enviroment_df = enviroment_df.drop_duplicates(subset=['name']) #Just in case
enviroment_df = enviroment_df.reset_index(drop=True)


#df = pd.read_csv(cinput,sep='\t')
#genome_df = df[ df[ 'type' ].isin(['genome']) ]          #dataframe of genome (to get size)
#genome_df = genome_df.reset_index(drop=True)
#genome_df  = df[ df[ 'type' ].isin(['gene']) ]            #dataframe of genes (to get orientations...)
#genome_df = genome_df.reset_index(drop=True)
#o_df  = df[ df[ 'type' ].isin(['NAP', 'RNAP', 'EXT', 'origin']) ]   #dataframe of objects (NAPs and bound RNAPs)

#Get sizes and number of objects
#------------------
nbp =  circuit_df['size'][0]       #number of bp in circuit
n_genome = len(genome_df['start']) #number of functional sites in genome (genes, operators, etc...)
N = len(o_df['pos'])               #number of objects bound to the DNA (NAPs, RNAPs,...)
n_env = len(enviroment_df['type']) #number of molecule types in the enviroment

#Check if structure is linear or circular
#------------------
if circuit_df['struc'][0] == 'circle':
    circular=True
elif circuit_df['struc'][0] == 'linear':
    circular=False
else:
    print( "Circuit's structure not recognized, assuming linear structure")
    circular=False

#Let's add a direction column for genome
#------------------
genome_df['direction'] = 0 #genes
for i in range(n_genome):
    if genome_df.iloc[i]['type'] != "gene": #Only genes have direction!
        continue
    if genome_df.iloc[i]['start'] < genome_df.iloc[i]['end']:
        genome_df.at[i,'direction'] = 1
    if genome_df.iloc[i]['start'] > genome_df.iloc[i]['end']:
        genome_df.at[i,'direction'] = -1

#Here, instead of this, maybe I should check that it is actually attacthed to a site in the genome.
#------------------
#And since we are here, let's add START and END to o_df
o_df['start'] = 0 #NANs
o_df['end'] = 0
for i in range(N):
    if not ( genome_df['name'].eq( o_df.iloc[i]['bound'] ) ).any():
        if not o_df.iloc['bound'] == 'DNA': #Enzymes such as topo's could just bind anywhere in the DNA
                                           #so we forgive them and do not count them
            print( "An object is not related to a genome site.", o_df.iloc[i]['bound']  )

    b = o_df.iloc[i]['bound']
    mask = genome_df['name'] == b
    o_df.at[i,'start'] = genome_df[mask]['start']
    o_df.at[i,'end'] = genome_df[mask]['end']

#Let's print some info of the system
#----------------------------------------------------------------
if circular:
    print("Running {0} frames for circular system of {1} bp \n".format(frames, nbp))
else:
    print("Running {0} frames for linear system of {1} bp \n".format(frames, nbp))

#Sort models to use
#----------------------------------------------------------------

#For the genome
genome_df['fun'] = ''    #This one is the initiation function
genome_df['params'] = 0  #This o
genome_df['params'] = genome_df['params'].astype(object)
for i in range(n_genome):
    if genome_df.iloc[i]['model'] == 'none':
        continue
    if genome_df.iloc[i]['model'] == 'poisson':
        genome_df.at[i,'fun'] = sm.P_binding_Poisson
    elif genome_df.iloc[i]['model'] == 'maxmin':
        genome_df.at[i,'fun'] = sm.promoter_curve_opening_E_maxmin
    elif genome_df.iloc[i]['model'] == 'eff_ener':
        genome_df.at[i,'fun'] = sm.promoter_curve_opening_E 
    elif genome_df.iloc[i]['model'] == 'sam':
        genome_df.at[i,'fun'] = sm.promoter_curve_Meyer
    else:
        print("Sorry, couldn't recognise initiation model for genome site:", genome_df.iloc[i]['name'])

    #And the params
    if genome_df.iloc[i]['oparams'] != 'none':
        genome_df.at[i, 'params' ] = np.loadtxt( genome_df.iloc[i]['oparams'] )

#And now the enviroment
enviroment_df['fun'] = '' #This will be erased at the end
for i in range(n_env):
    if enviroment_df.iloc[i]['model'] == 'none':
        continue
    if enviroment_df.iloc[i]['model'] == 'topo_continium':
        enviroment_df.at[i,'fun'] = sm.topoisomerase_activity
    elif enviroment_df.iloc[i]['model'] == 'gyrase_continium':
        enviroment_df.at[i,'fun'] = sm.gyrase_activity
    else:
        print("Sorry, couldn't recognise initiation model for enviroment:", enviroment_df.iloc[i]['name'])

#    if o_df.iloc[i]['type'] == 'NAP':
#        o_df.at[i,'direction'] = 0
#    elif o_df.iloc[i]['type'] == 'RNAP':
#        if o_df.iloc[i]['start'] < o_df.iloc[i]['end']:
#            o_df.at[i,'direction'] = 1
#        elif o_df.iloc[i]['start'] > o_df.iloc[i]['end']:
#            o_df.at[i,'direction'] = -1
#        elif o_df.iloc[i]['start'] == o_df.iloc[i]['end']: #This is unlikely, but if RNAP is in the termination?
#                                                           #I'll just stop it and I'll fix it later...
#            print("An RNAP in its end? This is not wrong, but I'll stop it. Fix your code Victor! and sorry user")
#            sys.exit()


#---------------------------------------------------------------------------------------------------------------------
#PREPARE INITIAL SYSTEM
#---------------------------------------------------------------------------------------------------------------------

#TENGO QUE ANADIR OTRA PARTE DISTINTA, QUE TRATE EL INPUT SISTEM DEPENDIENDO SI ES CONTINUACION O NO.
#ESTO MI PROGRAMA LO PUEDE SABER DEPENDIENDO SI ESTAN LOS EXT EN LOS INPUTS.
#TENGO QUE VOLVER A REVISAR ESTA PARTE ENTONCES

#Let's add fake ends (to stablish boundaries for either linear or circular structures)
#----------------------------------------------------------------
if continuation:    #I need to fix this. It would be better if the output doesn't have EXT_L and EXT_R
    a = 'EXT_L' in o_df['name'].unique()
    b = 'EXT_R' in o_df['name'].unique()
    print('Resuming simulation')
    if not (a and b):
        print('There is something wrong with the continuation file')
        print('Bye bye')
        sys.exit()

else:        #If it is a new run
    if N > 0: #This can only happen if there are objects bound to the DNA
        o_df['superhelical'] =circuit_df.iloc[0]['superhelical'] #I added this on 17/08/2022 (I think this is 
                                                                # correct)

        if circular: #For circular DNA
            start00, startNN = mm.get_start_end_c( o_df.iloc[0], o_df.iloc[-1], nbp) 
            #s_0 = o_df.iloc[-1]['superhelical']           #and the superhelical density
            s_0 = circuit_df.iloc[0]['superhelical']
            s_N = 0 #there's no point in using it (we won't need it)
            #t_0 = o_df.iloc[-1]['twist']
            t_0 = circuit_df.iloc[0]['twist']
            t_N = 0 # we don't need it here

        else: #For linear DNA
            if N > 1:
                s_0 = o_df.iloc[1]['superhelical'] #just because 
            else: #if it's one?
                s_0 = 0
                t_0 = 0
            s_N = 0 #there's no point in using it (we won't need it)
            t_N = 0
            start00   = 1
            startNN   = nbp

    else: #If nothing is bound
        start00 = 1
        startNN = nbp #it is the same in this case for either linear or circular
        s_0 = circuit_df.iloc[0]['superhelical']
        s_N = 0
        t_0 = circuit_df.iloc[0]['twist']
        t_N = 0

    #Let's treat the boundaries of our system as objects.
    #----------------------------------------------------------------
    #Notice that we don't specify the end
    extra = pd.DataFrame().reindex_like(o_df) #An empty copy
    extra = extra.drop_duplicates(subset=['name']) #This will only give us one row.
    extra = extra.append(extra, ignore_index=True) #Now we have two rows
    extra['type'] = ''
    extra['name'] = ''
    extra['size'] = 0
    extra['direction'] = 0
    #Let's fill the data of this extra
    extra.at[0,'start']        = start00  #Left boundary
    extra.at[0,'pos']          = start00
    extra.at[0,'superhelical'] = s_0
    extra.at[0,'twist']        = t_0
    extra.at[0,'type']         = 'EXT'       
    extra.at[0,'name']         = 'EXT_L'        
    extra.at[1,'start']        = startNN #Boundary on the right
    extra.at[1,'pos']          = startNN
    extra.at[1,'superhelical'] = s_N
    extra.at[1,'twist'] =        t_N
    extra.at[1,'type']         = 'EXT'       
    extra.at[1,'name']         = 'EXT_R'        
    #WARNING!!!!
    #There could be a big mistake in case of linear structures that have a NAP in positions 1 or nbp

    o_df = o_df.append(extra, ignore_index=True) #append left and right ends
    o_df = o_df.sort_values(by=["start"]) #and sort
    o_df = o_df.reset_index(drop=True)

N = len(o_df['start'])                #and update the number of objects
#----------------------------------------------------------------

#Ok, so even if twist is included, we assume that we don't know if we
#are resuming a simulation or if it's a new one.
#----------------------------------------------------------
#print(o_df)
for i in range(N-1): #And calculate the true twist
    #twist = mm.calculate_twist( o_df.iloc[i], o_df.iloc[i+1] )
    #o_df.iloc[i]['twist'] = twist
#    print(i)
    o_df.at[i,'twist'] = mm.calculate_twist( o_df.iloc[i], o_df.iloc[i+1] )

#---------------------------------------------------------------------------------------------------------------------
#START
#---------------------------------------------------------------------------------------------------------------------
#The motion starts, each RNAP moves with constant velocity v0, inyecting gamma supercoiling ahead and behind.

#print(o_df)
#print("")
#sys.exit()

sfile = open(sigmafile, "w")
pfile = open(posfile, "w")
nfile = open(namefile, "w")

#Let's quickly write the first parameters
#--------------------------------------------------------------------------
pos_l  = [] #since we are here, let's record the positions
sup_l  = [] #supercoiling
name_l = [] #and names
for i in range(0,N-1): 
    pos_l.append( o_df.iloc[i]['pos'] ) #record position
    sup_l.append( o_df.iloc[i]['superhelical'] ) #record supercoiling
    name_l.append( o_df.iloc[i]['name'] ) #record name
pos_l  = '  '.join( str(n) for n in pos_l) + '\n'
sup_l  = '  '.join( str(n) for n in sup_l) + '\n'
name_l = '  '.join( str(n) for n in name_l) + '\n'

pfile.write(  pos_l )
sfile.write(  sup_l )
nfile.write( name_l )

#additional*****
globalsigma = np.zeros( (frames,2) ) #this will save the global supercoiling level 
                            #at each time
#--------------------------------------------------------------------------

#print(genome_df)
#print("")
#sys.exit()
#control variables to corraborate our model is working
bind_rate = np.zeros( n_genome ) #this measures the binding times (to measure rate)
nb    = np.zeros( n_genome ) #to calculate the number of bound proteins

output_array = np.zeros(( frames, n_genome, 5 )) #this counts 
                                                #(1) the number of RNAPs reading the genes at each frame, 
                                                #(2) a array with 0 & 1 which indicate the time when a
                                                #RNAP binds each gene (0=no binding, 1=bind), 
                                                #(3) 0 or 1 if elongation terminated 
                                                #(4) the supercoiling & 
                                                #(5) twist  of each promoter
aux_array = np.zeros( (n_genome,5) ) #This one is an auxiliar of the one above but it is an auxiliar to fill it
                                    #So the format is: number of RNAPs elongating, binding, 
                                    #supercoiling and twist at promoter
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#FRAME BY FRAME
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
print("Sum of twist")
print(o_df['twist'].sum() )
print(o_df)

mask=o_df['type'] == 'RNAP'
n_RNAPs = len(o_df[mask]['start']) #Number of RNAPs bound

#print(genome_df)
#sys.exit()
for k in range(frames): #Go by frame 

    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #STATISTICAL PART!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    for i in range(n_genome): #Let's go gene by gene
        #se calcula la probabilidad de que un RNAP lo bindee
        #se echa un volado para ver si se bindea el RNAP
        #pero antes se debe de checar si el gen este disponible, o si no hay 
        #RNAP que bloquee el paso

       #This is just to ensure that the binding happens once...
       # if k  != 30:
       #     break

        #For now, only genes!
        #-----------------------------------------------------------
        if genome_df.iloc[i]['type'] != 'gene': 
            continue
        #-----------------------------------------------------------

        #CHECK IF GENE IS AVAILABLE
        #-------------------------------------------------------------
        start_i = genome_df.iloc[i]['start']        #the start site
        end_i   = genome_df.iloc[i]['end']          #the end site
        d_i     = genome_df.iloc[i]['direction']    #the direction
        mask=o_df['start'] <= start_i       #this is to locate the domain that previously contained
                                            #the twist of the partitioned region
        idx = o_df[mask].index[-1]          #the index of the domain before the binding
        if d_i == 1:
            d = o_df.iloc[idx+1]['pos'] - start_i #distance between promoter and RNAP on the right
        elif d_i == -1:
            d = start_i - o_df.iloc[idx]['pos']  #distance between promoter and RNAP on the left
        else:
            print('Error in checking if gene is occupied')
            print('There cannot be a gene without direction')
            sys.exit()


        #Let's record the levels of supercoiling and twist before binding...
        aux_array[i,3] = o_df.iloc[idx]['superhelical']
        aux_array[i,4] = o_df.iloc[idx]['twist']
        aux_array[i,1] = 0  #Let's reset the initiation to 0, if the RNAP binds then it'll be changed to 1
        aux_array[i,2] = 0  #And reset the termination to 0 as well

        if d <= o_df.iloc[idx+1]['size']:
            break #If the size of the new RNAP does not fit in d, then we skip the calculation as
                  #the RNAP cannot bind

        #WILL THE RNAP BIND?
        #-------------------------------------------------------------
        #According to the rate, calculate the probability of binding

        #(Probably here, it'll be way more useful/easy to use classes)
        prob = sm.P_binding( genome_df.iloc[i], o_df.iloc[idx]['superhelical'])
        
        #Now, let's decide if the RNAP will bind
        uran = np.random.uniform() #we need a random number
        
        #DECIDE
        #-----------------------------------------
        if uran <= prob:           #and decide
            #print("RNAP bound \n")

            #update control variables
            bind_rate[i] += k*dt-bind_rate[i]
            nb[i]    += 1

            aux_array[i,0] += 1 #Let's count the RNAP that just bound
            aux_array[i,1] = 1  #And let's register the time of initiation
            #print( i, aux_array[i,0] )

            #CREO QUE ANTES DE ACTUALIZAR ESTO, SE DEBE DE PRIMERO ESTABLECER LOS LENGTHS Y TODO ES
            #OSEA, CREO QUE ES ANADIRLA, ESTABLECER OTRA VEZ LOS LENGTHS, PARA EL CASO CIRCULAR,
            #Y LUEGO DECIFRAR EL TWIST
            #TAL VEZ ES MEJOR CON UN ATRIBUTO

            #Add RNAP
            #--------------------------------------------------------
            #These quantities define the start and end site of the new
            #RNAP + its direction

            #For now, only RNAPs can bind but in the future, I'll need to find a way to
            #figure out how the other enzymes will recognize the sites and bind
            #-----------------------------------------------------------
            mask = enviroment_df['site'] == genome_df.iloc[i]['type']
            size = enviroment_df[mask]['size'].to_string(index=False)
            stype = enviroment_df[mask]['type'].to_string(index=False) 
            sname = enviroment_df[mask]['name'].to_string(index=False)

            #Add the new object
            #-----------------------------------------------------------
            extra = pd.DataFrame().reindex_like(extra) #An empty copy
            extra = extra.drop_duplicates(subset=['name'], ignore_index=True) #This will only give us one row.
            extra['type'] = ''
            extra['name'] = ''
            extra['bound'] = ''
 
            extra.at[0, 'type'] = stype
            extra.at[0, 'name'] = sname
            extra.at[0, 'bound'] = genome_df.iloc[i]['name']
            extra.at[0, 'size'] = size
            extra.at[0, 'direction'] = d_i
            extra.at[0,'start'] = start_i
            extra.at[0,'end'] = end_i
            extra.at[0,'pos'] = start_i #it's position is where it starts (obviusly)
            extra.at[0,'twist'] = 0
            extra.at[0,'superhelical'] = 0

            T_twist     = o_df.iloc[idx]['twist'] #The twist and length (previous binding)
            T_length    = mm.calculate_length( o_df.iloc[idx], o_df.iloc[idx+1] )
            o_df = o_df.append(extra, ignore_index=True) #And add the new domain
            o_df = o_df.sort_values(by=["start"]) #sort
            o_df = o_df.reset_index(drop=True)
            N = len(o_df['start'])                #and update the number of objects

            #We need to update the positions of the fake boundaries in circular DNA
            #--------------------------------------------------------------------------
            if circular: #For circular DNA

                start00, startNN = mm.get_start_end_c( o_df.iloc[1], o_df.iloc[N-2], nbp) 
                o_df.at[0,'start']  = start00
                o_df.at[N-1,'start'] = startNN
                o_df.at[0,'pos']  = start00
                o_df.at[N-1,'pos'] = startNN

            #We are still missing the supercoiling density and the excess of twist...
            #We need to partition the twist so it is conserved...
            #--------------------------------------------------------
            length_l    = mm.calculate_length( o_df.iloc[idx], extra.iloc[0])
            length_r    = mm.calculate_length( extra.iloc[0], o_df.iloc[idx+2] ) #it is +2 because the
                                                                                 #RNAP has been inserted
            #So length_l, will be the new domain on the left and length_r the one on the right
            #now to calculate the new twists
            twist_l = T_twist*(  ( length_l + 0.5*extra.iloc[0]['size'] )/T_length  )
            twist_r = T_twist*(  ( length_r + 0.5*extra.iloc[0]['size'] )/T_length  )
            
            #update twists
            #------------CIRCULAR DNA--------------------
            if circular:  

                #There is no other domains besides the newly bound protein.
                if o_df.iloc[idx]['name'] == 'EXT_L' and o_df.iloc[idx+2]['name'] == 'EXT_R':
                    twist_i = T_twist
                    #In this case, the twist of EXT_L and twist_i remain the same
                    #becauwe it is a circular DNA with only one RNAP (no NAPs)

                #There is one EXT at the left
                elif o_df.iloc[idx]['name'] == 'EXT_L' and o_df.iloc[idx+2]['name'] != 'EXT_R':
                    o_df.at[idx, 'twist'] = twist_l 
                    o_df.at[N-2, 'twist'] = twist_l #The last one is equal to EXT_L
                    twist_i = twist_r

                #There is one EXT at the right
                elif o_df.iloc[idx]['name'] != 'EXT_L' and o_df.iloc[idx+2]['name'] == 'EXT_R':
                    o_df.at[idx, 'twist'] = twist_l
                    o_df.at[0, 'twist']   = twist_r 
                    o_df.at[N-2, 'twist'] = twist_r #The last ones are connected
                    twist_i = twist_r

                #In any other case
                else:
                    o_df.at[idx,'twist'] = twist_l         #update the one on the left
                    twist_i = twist_r

            #------------LINEAR DNA--------------------
            else:       
                o_df.at[idx,'twist'] = twist_l         #update the one on the left
                twist_i = twist_r

            o_df.at[idx+1, 'twist'] = twist_i #Update the twist of the new RNAP

            #--------------------------------------------------------

            for i in range(N-1): #And calculate the true supercoiling
                o_df.at[i,'superhelical'] = mm.calculate_supercoiling( o_df.iloc[i], o_df.iloc[i+1] )

            mask=o_df['type'] == 'RNAP'
            n_RNAPs = len(o_df[mask]['start']) #Number of RNAPs bound

    #--------------------------------------------------------
    #--------------------------------------------------------


    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #MECHNICAL PART!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------

    #THE MOTION AND TWIST INJECTION
    #--------------------------------------------------------------------------
    #if N_RNAPs > 0: #There is no point in doing this loop if there are no RNAPs bound
    for i in range(1,N-1):    #Let's go through each object, updating its twist and
                               #We ignore the last fake NAP as it's overtwist won't change
        #Notice that if N<=2, this loop never happens because there are not RNAPs
        #bound to the DNA

        #MOVING AN TWISTING
        #--------------------------------------------------------------------------
        #Is this an RNAP? (If not we can skip the loop)
        if o_df.iloc[i]['direction'] != 0:

            #Ok, it is a RNAP, so the motion starts and...
            #Object i MOVES!
            #-----------------
            o_df.at[i,'pos'] = mm.RNAP_motion( o_df.iloc[i] )

            #SUPERCOILING is injected and modifies TWIST
            #-----------------
            #TWIST injected on the left and right:
            twist_left, twist_right = mm.twist_injected( o_df.iloc[i] )

            o_df.at[i,'twist']  += twist_right #The one on the right always affects itself

            if   circular and i == 1:  #In case if we affect the boundary on the left
                o_df.at[N-2,'twist'] += twist_left
            #elif circular and i == N-2: #Boundary on the right
            #    o_df.at[i-1,'twist'] += twist_left
                #o_df.at[1,'twist']   += twist_right
            else:                       #If linear or not at the ends in circular
                o_df.at[i-1,'twist'] += twist_left

        #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------

    #IF IT REACHES THE END?
    #--------------------------------------------------------------------------
    #for i in range(1,N-1):    #We need to go through all objects again    
    i = 0
    for j in range(1,N-1):    #We need to go through all objects again    
        i = i+1

        #So, when a RNAP is removed, the domain on the left absorves the twist...

        #condition for transcription in >>>>>>>>>>>>> right direction or
        #condition for transcription in <<<<<<<<<<<<< left  direction
        if ( o_df.iloc[i]['direction'] == 1  and o_df.iloc[i]['end'] - o_df.iloc[i]['pos'] <= 0 ) or \
           ( o_df.iloc[i]['direction'] == -1 and o_df.iloc[i]['end'] - o_df.iloc[i]['pos'] >= 0 ):

            #print("frame", k)
            #print(o_df)
            #print(o_df['twist'].sum() )

            #update twists -  because i t is absorved
            #------------CIRCULAR DNA--------------------
            if circular:  

                #There is no other domains besides the newly bound protein.
                if o_df.iloc[i-1]['name'] == 'EXT_L' and o_df.iloc[i+1]['name'] == 'EXT_R':
                    o_df.at[i, 'twist'] = o_df.iloc[i]['twist']
                    o_df.at[0,'twist'] = o_df.at[N-2,'twist']  #I added this on 15/08/2022

                #There is one EXT at the left
                elif o_df.iloc[i-1]['name'] == 'EXT_L' and o_df.iloc[i+1]['name'] != 'EXT_R':
                    o_df.at[N-2,'twist'] += o_df.at[i,'twist'] #twist on the other side is absorbed
                    o_df.at[0,'twist'] = o_df.at[N-2,'twist']  #I added this on 15/08/2022

                #In any other case
                else:
                    o_df.at[i-1,'twist'] += o_df.at[i,'twist'] #twist is absorved by the domain on the left

            #------------LINEAR DNA--------------------
            else:       
                o_df.at[i-1,'twist'] += o_df.at[i,'twist'] #twist is absorved by the domain on the left
            #------------------------------------------

            #Before dropping, let's find which gene it just transcribed
            mask = genome_df['end'] == o_df.iloc[i]['end']
            a = genome_df[mask].index
            aux_array[a,0] -=1                #And let's remove the unbound RNAP from the count
            aux_array[a,2] = 1  #And let's register the time of termination

            #Drop object and update number of objects/RNAPs
            o_df = o_df.drop([i])                      #And we drop the object
            o_df = o_df.reset_index(drop=True)         #And reset indixees
            N = len(o_df['start'])                     #and update the number of objects

            mask=o_df['type'] == 'RNAP'
            n_RNAPs = len(o_df[mask]['start']) #Let's count the number of RNAPs bound
            i = i-1 #so we keep tracking the protein because if now N, changed, but it didn't affect
                    #the loop

    #--------------------------------------------------------------------------

    #We need to update the positions of the fake boundaries in circular DNA
    #--------------------------------------------------------------------------
    if circular: #For circular DNA
        #if n_RNAPs > 0:
        #if N > 0:
        if N > 2:
            start00, startNN = mm.get_start_end_c( o_df.iloc[1], o_df.iloc[N-2], nbp) 
            o_df.at[0,'start']  = start00
            o_df.at[N-1,'start'] = startNN
            o_df.at[0,'pos']  = start00
            o_df.at[N-1,'pos'] = startNN
        else:
            o_df.at[0,'start']  = 1
            o_df.at[N-1,'start'] = nbp
            o_df.at[0,'pos']  = 1
            o_df.at[N-1,'pos'] = nbp
    #--------------------------------------------------------------------------

    #UPDATE SUPERCOILING
    #--------------------------------------------------------------------------
    for i in range(0,N-1):  
        o_df.at[i,'superhelical'] = mm.calculate_supercoiling( o_df.iloc[i], o_df.iloc[i+1] )
        #Now it'll be ready for the next iteration
    #print(o_df)
    #sys.exit()
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #TOPOISOMERASE I AND GYRASE ACTIVITY
    #--------------------------------------------------------------------------
    #--------------------------------------------------------------------------
    #print(k)
    for i in range(0,N-1):
        #these two quantities are the amount of supercoiling each enzyme removes
        topo_sigma = sm.topoisomerase_activity( o_df.iloc[i]['superhelical'] )
        gyra_sigma = sm.gyrase_activity(        o_df.iloc[i]['superhelical'] )
        o_df.at[i,'superhelical'] += topo_sigma + gyra_sigma #and add the contributions
        o_df.at[i,'twist'] = mm.calculate_twist( o_df.iloc[i], o_df.iloc[i+1] ) #and update twice since we're here
        #note that it is correct to perform this twist update as the only parameter that we need from the
        #next 'object' it is its position
    #--------------------------------------------------------------------------

    #And update the twist and supercoiling density of left boundary
    #--------------------------------------------------------------------------
    if circular:
        o_df.at[0,'twist']  = o_df.iloc[N-2]['twist']
        o_df.at[0,'superhelical']  = o_df.iloc[N-2]['superhelical']
    #--------------------------------------------------------------------------

    #Let's quickly write the first parameters
    #--------------------------------------------------------------------------
    pos_l  = [] #since we are here, let's record the positions
    sup_l  = [] #supercoiling
    name_l = [] #and names
    for i in range(0,N-1): 
        pos_l.append( o_df.iloc[i]['start'] ) #record position
        sup_l.append( o_df.iloc[i]['superhelical'] ) #record supercoiling
        name_l.append( o_df.iloc[i]['name'] ) #record name
    pos_l  = '  '.join( str(n) for n in pos_l) + '\n'
    sup_l  = '  '.join( str(n) for n in sup_l) + '\n'
    name_l = '  '.join( str(n) for n in name_l) + '\n'

    pfile.write(  pos_l )
    sfile.write(  sup_l )
    nfile.write( name_l )
    #--------------------------------------------------------------------------

    #Additional*****************
    #Let's save the global superhelical density and time (to observe the evolution)
    #--------------------------------------------------------------------------
    if circular:
        if N > 2:
            globalsigma[k,0] = (o_df['twist'].sum() - o_df.iloc[0]['twist'])/(w0*nbp)
        else:
            globalsigma[k,0] = o_df['twist'].sum()/(w0*nbp)
    else: #linear
        globalsigma[k,0] = o_df['twist'].sum()/(w0*nbp)
    globalsigma[k,1] = k*dt
    #--------------------------------------------------------------------------

    #And COUNT RNAPs transcribing, initition time, termination time &
    #twist/supercoiling at promoter
    #--------------------------------------------------------------------------
    output_array[k,:,:] = aux_array

#--------------------------------------------------------------------------

#CLOSE THE ANIMATION FILES
#--------------------------------------------------------------------------
pfile.close()
sfile.close()
nfile.close()

#and the additional global sigma
np.savetxt('global_sigma.txt', globalsigma)

#Let's save the output array
np.save('output.npy', output_array)

#OUTPUT DATAFRAMES
#--------------------------------------------------------------------------
#Let's first drop the directions (they're in the way as they can be infered...)
o_df = o_df.drop(columns='start')
o_df = o_df.drop(columns='end')

#o_df = o_df.drop(columns='direction')
genome_df = genome_df.drop(columns='direction')
genome_df = genome_df.drop( columns='fun'  )
genome_df = genome_df.drop( columns='params'  )

enviroment_df = enviroment_df.drop( columns='fun'  )

#Circuit will now contain the average supercoiling and twist.
print(N)
if circular:
    if N > 2:
        circuit_df.at[0, 'superhelical'] = (o_df['twist'].sum() - o_df.iloc[0]['twist'])/(w0*nbp)
        circuit_df.at[0, 'twist'] = (o_df['twist'].sum() - o_df.iloc[0]['twist'])
    else:
        circuit_df.at[0, 'superhelical'] = o_df['twist'].sum()/(w0*nbp)
        circuit_df.at[0, 'twist'] = o_df['twist'].sum()
else: #linear
    circuit_df.at[0, 'superhelical'] = o_df['twist'].sum()/(w0*nbp)
    circuit_df.at[0, 'twist'] = o_df['twist'].sum()

#Print objects output -  in the future the enviroment may change as well...
o_df.to_csv('objects-'+coutput+'.csv', index=False, sep='\t')

#PRINT INFORMATION
#--------------------------------------------------------------------------
print("")
print("")
print("----------------------------------")
print("Simulation finished")
print("----------------------------------")
print("Average rates:")
print(bind_rate)
print(bind_rate/nb)
print("")
print("Proteins bound")
print(nb)
print("Sum of twist")
print(o_df['twist'].sum() )
print("Total superhelical density")
print(circuit_df.iloc[0]['superhelical'] )



