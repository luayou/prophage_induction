#!/bin/bash
#vs2 2nd round

#SBATCH --job-name=vs2_2ndround
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

# === PATHS AND FILES ===

main=<PATH_TO_MAIN_FOLDER>
phage_genomes=<PATH_TO_SELECTED_GENOMES_AFTER_1ST_ROUND_CLEANING_FASTA>
out_folder=<PATH_TO_OUTFOLDER>

mkdir ${out_folder}
cd ${out_folder}

# === VirSorter2 parameters ===

module purge
module load <ANACONDA_MODULE>
source activate <PATH_TO_VS2_CONDA_ENV>

min_length=1000
groups="dsDNAphage,ssDNA"
threads=<NUM_CORES>

virsorter run \
	--seqname-suffix-off \
	--viral-gene-enrich-off \
	--prep-for-dramv \
	-w ${out_folder}/vs2_2ndround \
	-i ${phage_genomes} \
	--min-length ${min_length} \
	--include-groups ${groups} \
	-j ${threads} \
	all

echo "Done virsorter"

#=== DRAM ===

module purge
module load <ANACONDA_MODULE>
source activate <PATH_TO_DRAM_CONDA_ENV>

input=${out_folder}/vs2_2ndround/for-dramv

DRAM-v.py annotate -i ${input}/final-viral-combined-for-dramv.fa -v ${input}/viral-affi-contigs-for-dramv.tab -o ${out_folder}/dram --threads 24 --min_contig_size 1000

DRAM-v.py distill -i ${out_folder}/dram/annotations.tsv -o ${out_folder}/dramv-distill


echo "done dram"


