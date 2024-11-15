#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-15

#Create array with the sample IDs
FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")
XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

#Move to directory containing the sam files for each sample
SAM_DIR="/data/users/lland/rna_seq/ref_genome/mapped_reads"
cd ${SAM_DIR}

#Run samtools view command to convert each sam file to a bam file, its binary equivalent
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools view -hbS ${XX}_mappedReads.sam > ${XX}_mappedReads.bam