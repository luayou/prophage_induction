#!/bin/bash
#Trimo phred 20 from raw reads

#SBATCH --job-name=trimmo
#SBATCH -p <PARTITION>
#SBATCH -N 1
#SBATCH -n <NUM_CORES>
#SBATCH --mem=<MEMORY>
#SBATCH --time=<TIME_LIMIT>
#SBATCH -o <PATH_TO_LOG_FILE>

SECONDS=0

# === FILES AND PATHS ===

reads=<PATH_TO_RAW_DATA_READS_FOLDER>
output=<PATH_TO_OUTPUT_FOLDER>


#=== GETTING UNIQUE PREFIX OF PE READS ===

cd $reads

ls *R1_001.fastq.gz > bases.txt


#=== TRIMMOMATIC ===

mkdir ${output}
cd ${output}

while read line; do
       java -Xmx32G -jar /usr/local/trimmomatic/0.38/trimmomatic-0.38.jar PE -threads 24 -phred33 -trimlog "${line}.log" -basein $reads/$line -baseout $line SLIDINGWINDOW:4:20 CROP:57 HEADCROP:11 MINLEN:30  ;
done < $reads/bases.txt


duration=$SECONDS
echo "$((duration / 3600)) hours
