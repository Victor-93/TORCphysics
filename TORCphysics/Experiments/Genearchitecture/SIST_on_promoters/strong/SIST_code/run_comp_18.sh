#!/bin/bash

#30/01/2024

#Description
#------------------------------------------------------------------------------------------------------
# The idea is to run SIST for a range of superhelicities.

#Initial conditions
#------------------------------------------------------------------------------------------------------
mydir=dir_18
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
perl master.pl -f $inputfasta -a M -o sidd_0.dat -b  -p -r -s -0.18
perl master.pl -f $inputfasta -a M -o sidd_1.dat -b  -p -r -s -0.181
perl master.pl -f $inputfasta -a M -o sidd_2.dat -b  -p -r -s -0.182
perl master.pl -f $inputfasta -a M -o sidd_3.dat -b  -p -r -s -0.183
perl master.pl -f $inputfasta -a M -o sidd_4.dat -b  -p -r -s -0.184
perl master.pl -f $inputfasta -a M -o sidd_5.dat -b  -p -r -s -0.185
perl master.pl -f $inputfasta -a M -o sidd_6.dat -b  -p -r -s -0.186
perl master.pl -f $inputfasta -a M -o sidd_7.dat -b  -p -r -s -0.187
perl master.pl -f $inputfasta -a M -o sidd_8.dat -b  -p -r -s -0.188
perl master.pl -f $inputfasta -a M -o sidd_9.dat -b  -p -r -s -0.189

rm input* one_line*

