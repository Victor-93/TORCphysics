#!/bin/bash

#08/05/2022

#Description
#------------------------------------------------------------------------------------------------------
#I need to prepare the data for plotting

#Initial conditions
#------------------------------------------------------------------------------------------------------
path=../Benham_code/
sigma0=-200 #Initial supercoiling density
sigmaf=0  #FINAL
ds=1000
nbp=214

#Begin process
#------------------------------------------------------------------------------------------------------

START=$sigma0
END=$sigmaf
s=0

cp $path*.txt .
for ((i=$START; i<=$END; i++))
do
        sigma=`echo "scale=3; $i / $ds" | bc`
	file=sidd_$s.dat
	cp $path$file .

	#Numbe rof lines in the file
	lines=($(wc $file))
	d=$((lines-nbp+1))  #We need the last 214 lines....
	                    #I think I could've also done this with cat or tail?

	#I rename the file
	sed -n "$d,$"p $file > test
	awk '{print $1, $3}' test > P_$s.txt
	awk '{print $1, $4}' test > G_$s.txt
	awk '{print $2}' test > sequence.txt
	s=$((s+1))
	#echo $path$file
done
rm *dat

