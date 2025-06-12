#!/bin/bash

#30/01/2024

#Description
#------------------------------------------------------------------------------------------------------
# I will extract the conditionals and put them in files. 

#Initial conditions
#------------------------------------------------------------------------------------------------------
path=../SIST_code/
nbp=552
sigma0=0 #Initial superhelical density
dsigma=-1
ds=1000

# For the directories
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
		
		# Cope file here.
        	sigma=`echo "scale=3; $sigma0 / $ds" | bc`
		file=sidd_$i.dat
		#cp $path$file .


		echo $s $sigma
		
		echo $sigma > val1

        	# P_M
        	#--------------------------------------------------------------
        	awk '$1 == "Transition" {print $4}' $file > val2
        	if [ "$s" -eq 0 ]; then
                	paste val1 val2 > P_M.txt
        	else
                	paste val1 val2 >> P_M.txt
        	fi

                s=$((s+1))
		sigma0=$((sigma0 - 1))
        done
done
rm *dat
