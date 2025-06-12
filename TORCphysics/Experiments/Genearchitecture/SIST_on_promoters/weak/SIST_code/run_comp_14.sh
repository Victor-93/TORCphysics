#!/bin/bash

#30/01/2024

#Description
#------------------------------------------------------------------------------------------------------
# The idea is to run SIST for a range of superhelicities.

#Initial conditions
#------------------------------------------------------------------------------------------------------
mydir=dir_14
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
perl master.pl -f $inputfasta -a M -o sidd_0.dat -b  -p -r -s -0.14
perl master.pl -f $inputfasta -a M -o sidd_1.dat -b  -p -r -s -0.141
perl master.pl -f $inputfasta -a M -o sidd_2.dat -b  -p -r -s -0.142
perl master.pl -f $inputfasta -a M -o sidd_3.dat -b  -p -r -s -0.143
perl master.pl -f $inputfasta -a M -o sidd_4.dat -b  -p -r -s -0.144
perl master.pl -f $inputfasta -a M -o sidd_5.dat -b  -p -r -s -0.145
perl master.pl -f $inputfasta -a M -o sidd_6.dat -b  -p -r -s -0.146
perl master.pl -f $inputfasta -a M -o sidd_7.dat -b  -p -r -s -0.147
perl master.pl -f $inputfasta -a M -o sidd_8.dat -b  -p -r -s -0.148
perl master.pl -f $inputfasta -a M -o sidd_9.dat -b  -p -r -s -0.149

rm input* one_line*

