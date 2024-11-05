#!/usr/bin/env bash
#do the sbatch commands

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=800M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8

#change this line to the path that contains the modules for fastaQC here
XX="SRR7821919"
module load FastQC/0.11.9-Java-11

mkdir /data/users/${USER}/rna_seq/fastQC_${XX}
cd /data/users/${USER}/rna_seq/fastQC_${XX}

#Create symbolic links to the content in /data/courses/rnaseq_course/toxoplasma_de/reads/SRR7821919_2.fastq.gz
ln -s /data/courses/rnaseq_course/toxoplasma_de/reads/${XX}_1.fastq.gz ${XX}_1.fastq.gz
ln -s /data/courses/rnaseq_course/toxoplasma_de/reads/${XX}_2.fastq.gz ${XX}_2.fastq.gz


#Create symbolic links to the content in /data/courses/rnaseq_course/toxoplasma_de/reads/SRR7821919_2.fastq.gz
#count how many of them there are; we want to run the scripts as many times as there are 

fastqc -t 2 ${XX}_*.fastq.gz
