#!/bin/bash
#mafft

#SBATCH --job-name=mafft_linsi
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

#=== FILES AND PATHS ===

main=<PATH_TO_MAIN_FOLDER>
complete_ci_proteins=<PATH_TO_COMPLETE_CI_LIKE_REP_FASTA>
output_folder=<PATH_TO_OUTFOLDER>
out_file=<FILENAME_OUTPUT_FASTA>

#=== MAFFT ====

mkdir ${output_folder}

cd ${output_folder}

module purge
module load mafft/7.526

linsi ${complete_ci_proteins} > ${output_folder}/${out_file}

echo "linsi DONE"
