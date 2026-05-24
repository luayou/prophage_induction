#!/bin/bash
# VirSorter2 - 1st round

#SBATCH --job-name=vs2_1stround
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

module purge
module load <ANACONDA_MODULE>
source activate <PATH_TO_CONDA_ENV>

# === Paths ===
output_folder=<PATH_TO_OUTPUT_FOLDER>
bacterial_genomes=<PATH_TO_BACTERIAL_GENOMES_DIR>
genomes_file=<PATH_TO_GENOMES_LIST_FILE>

# === VirSorter2 parameters ===
min_length=5000
min_score=0.5
groups="dsDNAphage,ssDNA"
threads=<NUM_CORES>

mkdir ${output_folder}
cd ${bacterial_genomes}

while read i; do
    genome=${i%.*}
    virsorter run \
        -w ${output_folder}/${genome} \
        -i ${i} \
        --keep-original-seq \
        --min-length ${min_length} \
        --min-score ${min_score} \
        --include-groups ${groups} \
        -j ${threads} \
        all  # Sullivan protocol: https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v2?step=3
    echo "${i} done"
done < ${genomes_file}

echo "All jobs done"
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
-- INSERT --
#!/bin/bash
#1st vs2

#SBATCH --job-name=vs2_1stround
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem=32gb
#SBATCH --time=160:00:00
#SBATCH -o PATH



module purge
module load anaconda/5.0.1-Python2.7-gcc5

source activate /scratch/xa90/lave/software/conda_envs/vs2


output_folder=~/2.vs2_first_round/remaining_106_sof_isolates
bacterial_genomes=~/1.bacterial_genomes/

genomes_file=~/Sof_remaining_106genomes.txt

mkdir ${output_folder}

cd ${bacterial_genomes}

while read i ;do

        genome=${i%.*}
        virsorter run -w ${output_folder}/${genome} -i ${i} --keep-original-seq --min-length 5000 --min-score 0.5 --include-groups dsDNAphage,ssDNA -j 24 all  #following Sullivan protocol https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v2?step=3

        echo ${i} done

done < ${genomes_file}


echo "done"
~
~
~
~
~
~
~
~
~
~
~
~
~
-- INSERT --

~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
~
-- INSERT --
#!/bin/bash
#1st vs2

#SBATCH --job-name=vs2_1stround
#SBATCH -p comp
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem=32gb
#SBATCH --time=160:00:00
#SBATCH -o PATH



module purge
module load anaconda/5.0.1-Python2.7-gcc5

source activate /scratch/xa90/lave/software/conda_envs/vs2


output_folder=~/2.vs2_first_round/remaining_106_sof_isolates
bacterial_genomes=~/1.bacterial_genomes/

genomes_file=~/Sof_remaining_106genomes.txt

mkdir ${output_folder}

cd ${bacterial_genomes}

while read i ;do

        genome=${i%.*}
        virsorter run -w ${output_folder}/${genome} -i ${i} --keep-original-seq --min-length 5000 --min-score 0.5 --include-groups dsDNAphage,ssDNA -j 24 all  #following Sullivan protocol https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v2?step=3

        echo ${i} done

done < ${genomes_file}


echo "done"
~
~
~
~
~
~
~
~
~
~
~
~
~
-- INSERT --
