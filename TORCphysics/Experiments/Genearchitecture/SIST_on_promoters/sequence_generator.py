import numpy as np
import sys
import random

#---------------------------------------------------------------------------------------------
#DESCRIPTION
# I want to create fake sequences so I can infer the promoter opening curves
# of the weak, medium and strong promoters.

#Modify the process depending the case
#---------------------------------------------------------------------------------------------

#FUNCTIONS
#---------------------------------------------------------------------------------------------
#Generates sequence of length 'length' and only using bases in the string 'bases'
#The sequence generated is random...
def random_dna_sequence(length,bases):
    return ''.join(random.choice(bases) for _ in range(length))

#Calculate the percentage of GC content in a sequence
def calculate_GC_content(sequence):
    total = len(sequence)
    count = 0
    for s in sequence:
        if s == 'G' or s == 'C':
            count=count+1
    content = count/total
    return 100*content

#TECHNICAL STUFF
#---------------------------------------------------------------------------------------------
n_fasta=79       #Number of characters per line in fasta file

#INPUTS
#---------------------------------------------------------------------------------------------
inter_length = 250

# Promoter sequences
weak_seq ='AAAAAGAGTATTGACTTCGCATCTTTTTGTACCTATAATGTGTGGATAGCGG'
medium_seq = 'TTGACATCAGGAAAATTTTTCTGCATAATTATTTCATATCAC'
strong_seq = 'TTGACATCGCATCTTTTTGTACCTATAATGTGTGGATAGAGT'

cases_labels = ['weak', 'medium', 'strong']

#LET'S BUILD THE ACTUAL SEQUENCE OF EACH PROMOTER REGION
#---------------------------------------------------------------------------------------------
flaking_sequence = random_dna_sequence( inter_length, 'G' ) # Creates sequence of inter_length length made of 'G's

weak_system = flaking_sequence + weak_seq + flaking_sequence
medium_system = flaking_sequence + medium_seq + flaking_sequence
strong_system = flaking_sequence + strong_seq + flaking_sequence

#And write it in fasta format
#---------------------------------------------------------------------------------------------
seq_list = [weak_system, medium_system, strong_system]
for j, case in enumerate(cases_labels):
    output_file = case+'.fasta'
    mysequence = seq_list[j]

    nf  = len(mysequence)/n_fasta       #number of lines in fasta file
    n = int(len(mysequence)/n_fasta)    #number of lines in fasta file

    f= open(output_file,'w')
    for i in range(n):
        f.write(mysequence[i*n_fasta:n_fasta*(i+1)] )
        f.write('\n')
    if nf>n: #Just in case we have more
        f.write(mysequence[ n*n_fasta:] )
        f.write('\n')
    f.close()


