#!/bin/bash

#SBATCH --job-name=iqtree
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

module purge
module load iqtree

#=== FILES AND PATHS ===

main=<PATH_TO_MAIN_FOLDER>
alignment=<PATH_TO_CLEAN_ALIGNMENT_FASTA>
iqtree_log=<FILENAME_IQTREE_LOG>

#=== IQTREE PARAMETERS ===

mkdir ${main}/conservative_bootstrap

cd ${main}/conservative_bootstrap

model_selection=MFP
bootstrap_n=1000
threads=<THREADS>

iqtree -s ${alignment} \
	 -m ${model_selection}\
	 -b ${bootstrap_n} \
     -nt ${threads} &> ${iqtree_log}


