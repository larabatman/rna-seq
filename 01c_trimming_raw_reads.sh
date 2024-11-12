#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-15

source /data/users/lfalquet/SBC07107_24/scripts/module.sh

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar

FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")
XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

module load fastp/0.23.4-GCC-10.3.0

WORK_DIR="/data/users/${USER}/rna_seq"
cd $WORK_DIR

mkdir -p ${WORK_DIR}/fastQC_results/trimmed
mkdir -p ${WORK_DIR}/fastQC_results/fastQC_reports

ln -s /data/courses/rnaseq_course/toxoplasma_de/reads/${XX}_1.fastq.gz ${XX}_1.fastq.gz
ln -s /data/courses/rnaseq_course/toxoplasma_de/reads/${XX}_2.fastq.gz ${XX}_2.fastq.gz

#Run Trimmomatic on paired-end (PE) sequences

#java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -threads 4 ${XX}_1.fastq.gz ${XX}_2.fastq.gz \
 # ${WORK_DIR}/fastQC_results/trimmed/${XX}_1trim.fastq.gz ${WORK_DIR}/fastQC_results/trimmed/${XX}_1unpaired.fastq.gz \
 # ${WORK_DIR}/fastQC_results/trimmed/${XX}_2trim.fastq.gz ${WORK_DIR}/fastQC_results/trimmed/${XX}_2unpaired.fastq.gz \
 # LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100

#Sadly, with this parameter, the quality of the reads was too low and everything was dropped 
#Input Read Pairs: 36220112 Both Surviving: 0 (0.00%) Forward Only Surviving: 0 (0.00%) Reverse Only Surviving: 0 (0.00%) Dropped: 36220112 (100.00%)
#Try relaxing the conditions:
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -phred33 -threads 4 ${XX}_1.fastq.gz ${XX}_2.fastq.gz \
  ${WORK_DIR}/fastQC_results/trimmed/${XX}_1trim.fastq.gz ${WORK_DIR}/fastQC_results/trimmed/${XX}_1unpaired.fastq.gz \
  ${WORK_DIR}/fastQC_results/trimmed/${XX}_2trim.fastq.gz ${WORK_DIR}/fastQC_results/trimmed/${XX}_2unpaired.fastq.gz \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:50

#Run the fastqc analysis on both reverse and forward sequences for each sample
fastqc -t 2 ${WORK_DIR}/fastQC_results/trimmed/${XX}_1trim.fastq.gz ${WORK_DIR}/fastQC_results/trimmed/${XX}_2trim.fastq.gz -o ${WORK_DIR}/fastQC_results/fastQC_reports
