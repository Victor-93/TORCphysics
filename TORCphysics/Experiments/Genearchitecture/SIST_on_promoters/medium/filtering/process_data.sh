#!/bin/bash

#03/06/2024

#Description
#------------------------------------------------------------------------------------------------------
#I need to prepare the data for plotting

#Initial conditions
#------------------------------------------------------------------------------------------------------
path=../SIST_code/
nbp=542

start_d=0
end_d=19
#Begin process
#------------------------------------------------------------------------------------------------------

s=0
for ((d=$start_d; d<=$end_d; d++))
do
	# Copy data files here.
	dir=dir_$d/
	cp $path$dir*dat . 

	for ((i=0; i<=9; i++))
	do
	
		#And extract info
		#--------------------------------------------
		file=sidd_$i.dat
		#Numbe rof lines in the file
		lines=($(wc $file))
		l=$((lines-nbp+1))  #We need the last nbp lines....
	                    	    #I think I could've also done this with cat or tail?

		#I rename the file
		sed -n "$l,$"p $file > test
		awk '{print $1, $3}' test > Prob_M_$s.txt
		awk '{print $1, $4}' test > G_M_$s.txt
		awk '{print $2}' test > sequence.txt
		s=$((s+1))
	done
done

rm *dat

