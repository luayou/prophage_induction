#!/bin/bash

#SBATCH --job-name=gtbd-tk
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>


# === FILES AND PATHS ===

bacterial_genomes=<PATH_TO_GENOMES_LIST_FILE>
output_dir=<PATH_TO_OUTPUT_FOLDER>

# === SETTING UP THE ENV ===

module purge
module load <ANACONDA_MODULE>
source activate <PATH_TO_CONDA_ENV>
conda env config vars set GTDBTK_DATA_PATH=<PATH_TO_GTDBTK_DB_LAST_RELEASE>

gtdbtk classify_wf --genome_dir ${bacterial_genomes} -x fna --out_dir ${output_dir} --cpus 24 --write_single_copy_genes --keep_intermediates --skip_ani_screen
