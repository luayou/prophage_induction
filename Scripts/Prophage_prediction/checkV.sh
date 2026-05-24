#!/bin/bash
#checkV

#SBATCH --job-name=checkV
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

# === Setting up the environment ===

module purge
module load <ANACONDA_MODULE>
source activate <PATH_TO_CHECKV_CONDA_ENV>
export CHECKVDB=<PATH_TO_CHECKV_DB>

# === User-defined paths ===
output_folder=<PATH_TO_OUTPUT_FOLDER>
vs2_output_genomes=<PATH_TO_VS2_OUTPUT_GENOMES_FASTA>

mkdir ${output_folder}
cd ${output_folder}

checkv end_to_end ${vs2_output_genomes} checkv_output -t 24
