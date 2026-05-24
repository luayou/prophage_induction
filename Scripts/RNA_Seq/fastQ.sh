#!/bin/bash
#fastqc_ 1x

#SBATCH --job-name=fastqc
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>


module load fastqc/0.12.1

# === FILES AND PATHS ===

reads_folder=<PATH_TO_READS_FOLDER>
output=<PATH_TO_OUTPUT_FOLDER>

mkdir ${output}

cd ${output}

suffix_all_reads=${SUFFIX}

fastqc ${SUFFIX} -o ${output} --noextract

echo "DONE :)"
