#!/bin/bash

#30/01/2024

#Description
#------------------------------------------------------------------------------------------------------
# The idea is to run SIST for a range of superhelicities.

#Initial conditions
#------------------------------------------------------------------------------------------------------
mydir=dir_12
inputfasta=input.fasta
#sigma0=-8 #Initial supercoiling density
#sigmaf=0  #FINAL

#Begin process
#------------------------------------------------------------------------------------------------------
# Create directory and copy SIST in the
mkdir $mydir
cd $mydir
mkdir trans_three
cd ..
cp trans_three/qsidd $mydir/trans_three
cp $inputfasta $mydir
cp irf308.linux.exe $mydir
cp IR_finder.pl $mydir
cp master.pl $mydir
cd $mydir

#And run qsidd
perl master.pl -f $inputfasta -a M -o sidd_0.dat -b  -p -r -s -0.12
perl master.pl -f $inputfasta -a M -o sidd_1.dat -b  -p -r -s -0.121
perl master.pl -f $inputfasta -a M -o sidd_2.dat -b  -p -r -s -0.122
perl master.pl -f $inputfasta -a M -o sidd_3.dat -b  -p -r -s -0.123
perl master.pl -f $inputfasta -a M -o sidd_4.dat -b  -p -r -s -0.124
perl master.pl -f $inputfasta -a M -o sidd_5.dat -b  -p -r -s -0.125
perl master.pl -f $inputfasta -a M -o sidd_6.dat -b  -p -r -s -0.126
perl master.pl -f $inputfasta -a M -o sidd_7.dat -b  -p -r -s -0.127
perl master.pl -f $inputfasta -a M -o sidd_8.dat -b  -p -r -s -0.128
perl master.pl -f $inputfasta -a M -o sidd_9.dat -b  -p -r -s -0.129

rm input* one_line*

