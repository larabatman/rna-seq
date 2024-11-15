#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=25000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-15

#Create an array with the sample IDs
FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")
XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

#Move to the directory containing the bam files
BAM_DIR="/data/users/lland/rna_seq/ref_genome/mapped_reads"
cd ${BAM_DIR}

#Run samtools sort command to sort the bam files, organizing the alignments based on their genomic position
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools sort -m 25000M -@ 4 -o ${XX}_sorted.bam -T temp ${XX}_mappedReads.bam