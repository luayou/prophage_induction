#!/bin/bash
#getting_15kb_flanked_regions

#SBATCH --job-name=flanked_regions
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

# === PATHS ===

main=<PATH_TO_MAIN_FOLDER>
phage_genomes=<PATH_TO_VS2_OUTPUT_FASTA>
trimmed_phages_with_flanked_regions=<PATH_TO_NEW_FASTA_FILE_WITH_FLANKED_REGIONS>

bed_phages_to_cut=<PATH_TO_BED_FILE_WITH_COORDINATES_TO_ADD_15KB_TRIMMED_PHAGES_BASED_ON_ORIGINAL_PHAGE_PREDICTION>

#=== GETTING FASTA BEDTOOLS ===

module purge

module load bedtools/2.31.1

bedtools getfasta -fi ${phage_genomes} -bed ${bed_phages_to_cut} -fo ${trimmed_phages_with_flanked_regions} #-name
