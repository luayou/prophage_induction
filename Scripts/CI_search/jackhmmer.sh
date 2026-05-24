#!/bin/bash

#SBATCH --job-name=jackhmmer
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

# === PATHS AND FILES ===

main=<PATH_TO_MAIN_FOLDER>
phages_proteins=<PATH_TO_AA_PHAGE_PROTEINS_FASTA>
outdir=<PATH_TO_OUTFOLDER>
domains_output=<PATH_TO_TSV_SUMMARY_FILE_TOP_DOMAINS>
sequences_output=PATH_TO_TSV_SUMMARY_FILE_TOP_SEQUENCES>
hmm_profiles_per_iter=<PATH_TO_PREFIX_HMM_PER_ITER>
alignments_per_iter=<PATH_TO_PREFIX_ALI_PER_ITER>

ci=<PATH_TO_LAMBDA_CI_UNIPROTID_P03034_FASTA>

# == JACKHMMER PARAMETERS ===

iter=7
incE=1e-10
incdomE=1e-10
threads=<THREADS>

mkdir ${outdir}

module purge
module load hmmer/3.4

jackhmmer -N ${iter} --incE ${incE} --incdomE ${1e-10} --cpu ${threads} --chkhmm ${outdir}/${hmm_profiles_per_iter} --chkali ${outdir}/{alignments_per_iter} --tblout ${outdir}/${sequences_output} --domtblout ${outdir}/${domains_output} ${ci} ${phages_proteins}
