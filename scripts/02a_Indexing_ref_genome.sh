#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8

cd /data/users/lland/rna_seq/ref_genome

#Unzip the primary assembly fasta file containing the reference genome
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

#Build the index for that reference genome
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif hisat2-build Mus_musculus.GRCm39.dna.primary_assembly.fa GRCm39_genome