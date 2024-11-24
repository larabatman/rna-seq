#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8

BAM_DIR="/data/users/lland/rna_seq/ref_genome/mapped_reads"
cd ${BAM_DIR}

#Take a lil bit of choromosme 8 fom one of the WT lung samples
#Take a piece of one of the WT lung samples: chromosome 8 
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools view -b SRR7821918_sorted.bam "8" > SRR7821918_subset_chr8.bam
#Index that fraction of bam file 
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif samtools index SRR7821918_subset_chr8.bam