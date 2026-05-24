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

# === FILES AND PATHS ===

main=<PATH_TO_GTBKTK_OUTPUT>
alignment=<PATH_TO_gtdbtk.bac120.user_msa.fasta.gz>

cd ${main}

iqtree -s ${main}/${alignment} -m MFP -bb 1000 -nt 12 &> iqtree_only_bacteroides_isolates.log
