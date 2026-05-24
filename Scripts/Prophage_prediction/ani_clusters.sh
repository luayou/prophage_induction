#!/bin/bash
#ANI_1st_round and checkV

#SBATCH --job-name=ANI_selecting
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

module purge

# === Paths ===
main=<PATH_TO_MAIN_FOLDER>
out_folder=<PATH_TO_OUTPUT_FOLDER>
trimmed_genomes=<PATH_TO_CHECKV_TRIMMED_PROPHAGES_FASTA>
original_genomes=<PATH_TO_VS2_OUTPUT_FASTA>
phage_genomes=<PATH_TO_NEW_PHAGE_CAT_FASTA>
blast_output=<BLAST_OUTPUT_FILENAME>
ani_output=<ANI_OUTPUT_FILENAME>
clusters_output=<ANI_CLUSTERS_OUTPUT_FILENAME>

scripts=<PATH_TO_CHECKV_ANI_SCRIPTS_FOLDER>

#=== Blast all-against-all ===

module load blast/2.15.0
mkdir ${out_folder}

cat ${trimmed_genomes} ${original_genomes} > ${phage_genomes}

makeblastdb -in ${phage_genomes} -dbtype nucl -out ${phage_genomes}

blastn -query ${phage_genomes} -db ${phage_genomes} -outfmt '6 std qlen slen' -max_target_seqs 10000 -out ${out_folder}/${blast_output} -num_threads 32

#=== ANI clustering ==== 

module purge
module load checkv/0.7.0

python ${scripts}/anicalc.py -i ${out_folder}/${blast_output} -o ${out_folder}/${ani_output}

#=== ANI and COV parameters ===

ani=95
min_tcov=30 # we used this low value for our initial step of cleaning 
min_qcov=0

python ${scripts}/aniclust.py --fna ${phage_genomes} --ani ${out_folder}/${ani_output} --out ${out_folder}/${clusters_output} --min_ani ${ani} --min_tcov ${min_tcov} --min_qcov ${min_qcov}
