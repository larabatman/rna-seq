#!/usr/bin/env bash


#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pibu_el8

#Load MultiQC to produce an aggregated report of our fastqc results (compare quality of reads across samples)
module load MultiQC/1.11-foss-2021a

FASTQC_RESULTS_DIR="/data/users/${USER}/rna_seq/fastQC_results"
MULTIQC_OUTPUT_DIR="/data/users/${USER}/rna_seq/multiQC_results"

#Create a directory to store the multiQC report
mkdir $MULTIQC_OUTPUT_DIR

#Run multiQC:
# -o: direct the output
# --interactive: to get an HTML report
# --title: to differentiate between cleaned and uncleaned reads
multiqc $FASTQC_RESULTS_DIR -o $MULTIQC_OUTPUT_DIR --interactive --title "MultiQC Report - Raw Data"
