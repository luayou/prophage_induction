#!/bin/bash
#CDHIT complete CI-like repressors

#SBATCH --job-name=cd-hit
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>


module load cdhit/4.8.1

#=== FILES AND PATHS ===

main=<PATH_TO_MAIN_FOLDER>
prefix_output=<PREFIX_OUTPUT>
ci_proteins=<PATH_TO_CI_COMPLETE_PROTEINS_FASTA>

#=== CD-HIT PARAMETERS ===

cd ${main}

module load cdhit/4.8.1

identity=1
word=5
mode=1 #slow mode
type=0 #local
cov_shorter=1
stop=0
threads=<THREADS>
mem=<MEMORY>

cd-hit -i ${ci_proteins} \
	-o ${main}/${prefix_output} \
	-c ${identity} \
	-n ${word} \
	-g ${mode} \
	-d ${stop} \
	-G ${type} \
	-aS ${cov_shorter} \
	-T ${threads} \
	-M ${mem}
