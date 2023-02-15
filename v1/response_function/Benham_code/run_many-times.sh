#!/bin/bash

#08/05/2022

#Description
#------------------------------------------------------------------------------------------------------
#The objective is to obtain the promoter response curve.
#To that end, I need to run multiple times the SIDD code.

#Initial conditions
#------------------------------------------------------------------------------------------------------
sigma0=-200 #Initial supercoiling density
sigmaf=0  #FINAL
input=..//Pleu500_14bp.fasta
ds=1000

#sigma0 and sigmaf will be divided by ds

#Begin process
#------------------------------------------------------------------------------------------------------

START=$sigma0
END=$sigmaf
s=0
for ((i=$START; i<=$END; i++))
do
	sigma=`echo "scale=3; $i / $ds" | bc`
	out=sidd_$s.dat
	echo "Doing sigma=$sigma"
	if [ $s -eq 0 ]
	then
		echo $s $sigma > info.txt
	else
		echo $s $sigma >> info.txt
	fi
	s=$((s+1))
	time perl master.pl -f $input -a M -o $out -b -p -r -s $sigma
	#echo "time perl master.pl -f $input -a M -o $out -b -p -r -s $sigma"
done

