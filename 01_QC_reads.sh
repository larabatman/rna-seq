#!/usr/bin/env bash
#do the sbatch commands

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-2

#change this line to the path that contains the modules for fastaQC here

FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")
FASTQ_DIR="/data/courses/rnaseq_course/toxoplasma_de/reads/"


#FASTQ_FILES=($(ls ${FASTQ_DIR}/*_1.fastq.gz | sed 's/_1\.fastq\.gz//g' ))
#FASTQ_FILES=($(ls ${FASTQ_DIR}*_1.fastq.gz | grep -oP '(\S+)_1\.fastq\.gz'))

XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

module load FastQC/0.11.9-Java-11

mkdir -p /data/users/${USER}/rna_seq/fastQC_${XX}
cd /data/users/${USER}/rna_seq/fastQC_${XX}

#Create symbolic links to the content in /data/courses/rnaseq_course/toxoplasma_de/reads/SRR7821919_2.fastq.gz
ln -s /data/courses/rnaseq_course/toxoplasma_de/reads/${XX}_1.fastq.gz ${XX}_1.fastq.gz
ln -s /data/courses/rnaseq_course/toxoplasma_de/reads/${XX}_2.fastq.gz ${XX}_2.fastq.gz

fastqc -t 2 ${XX}_1.fastq.gz ${XX}_2.fastq.gz
