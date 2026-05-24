#!/bin/bash
#Map trimmo reads to Human

#SBATCH --job-name=to_human
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

SECONDS=0

#=== FILES AND PATHS ===

human_genome=<PATH_TO_HUMAN_GENOME> #GCF_000001405.40_GRCh38.p14_genomic.fna.gz
trimmo_reads=<PATH_TO_TRIMMED_READS>
wo_human_reads=<PATH_TO_OUTPUT_FOLDER>
output=<PATH_TO_TEMP_BAM_FOLDER>

#=== MAPP READS TO HUMAN AND FILTER THEM OUT ===

mkdir ${output}
mkdir ${wo_human_reads}

module load bowtie2/2.5.4
module load samtools/1.19.3

threads=${THREADS}

bowtie2-build ${human_genome} ${human_genome}

cd ${trimmo_reads}

for i in $(ls *_1P.fastq.gz); do
        rv=${i%_*}_2P.fastq.gz;
        outFile=${i%_*}
        bowtie2 -p 12 -x ${human_genome} -1 ${i} -2 ${rv} 2>${output}/${outFile}.log | samtools sort -o ${output}/${outFile}.sorted.bam -@ ${threads}
        samtools collate -u -O ${output}/${outFile}.sorted.bam -@ 12 | samtools fastq -f12 -1 ${wo_human_reads}/${outFile}_trimmoNohuman_1P.fq -2 ${wo_human_reads}/${outFile}_trimmoNohuman_2P.fq -s ${wo_human_reads}/${outFile}_trimmoNohuman_U.fq -n -@ ${threads}
done

echo It took $SECONDS seconds

echo "DONE"
