#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=1-15

#Retrieve fastq files for each sample ID
FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")
XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

#Create directory variables 
INDEX_DIR="/data/users/lland/rna_seq/ref_genome"
FASTQ_DIR="/data/users/lland/rna_seq"
OUTPUT_DIR=${INDEX_DIR}/mapped_reads


cd ${INDEX_DIR}

#Create new directory for output files
mkdir -p mapped_reads

#Run HiSAT2 using the indexed reference genome, specifying RF strandedness, the fastq files for each sample, the sam output and the number of threads
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2 -x ${INDEX_DIR}/GRCm39_genome --rna-strandness RF -1 ${FASTQ_DIR}/${XX}_1.fastq.gz -2 ${FASTQ_DIR}/${XX}_2.fastq.gz -S ${OUTPUT_DIR}/${XX}_mappedReads.sam -p 4