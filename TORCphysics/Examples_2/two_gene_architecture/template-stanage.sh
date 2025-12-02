#!/bin/bash
#SBATCH --job-name=genearch_2gene
#SBATCH --time=00:10:00
# Request 16 gigabytes of real memory (mem) - We can request up to 120G
#SBATCH --mem=16G
#Number of CPUs
#SBATCH --cpus-per-task=51

# load the module for the program we want to run
module load Anaconda3/2022.05

# Load environment
source activate torc

#Run the program
python genearch_2gene_variations_parallel.py
