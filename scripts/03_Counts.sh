#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=02:00:00
#SBATCH --partition=pibu_el8
#SBATCH --array=0-15



#It appears that there is a small typo in the container pathway:
#/containers/apptainer/subread_2.0.1--hed695b0_0.sif

#This version of featureCounts does not have the following option
#featureCounts -p --countReadPairs -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam

#Create an array for job partitioning
FASTQ_FILES=("SRR7821918" "SRR7821919" "SRR7821920" "SRR7821921" "SRR7821922" "SRR7821937" "SRR7821938" "SRR7821939" "SRR7821949" "SRR7821950" "SRR7821951" "SRR7821952" "SRR7821953" "SRR7821968" "SRR7821969" "SRR7821970")
XX="${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}"

#Directory of the BAM files
BAM_DIR="/data/users/lland/rna_seq/ref_genome/mapped_reads"
cd ${BAM_DIR}

#Run  featureCounts with the following options:
# -T 4: four threads
# -p: paired-end
# -s 2: strandedness
# -t exon: -t to define the feature to count; here, we want to count exons
# -g gene_id: to group the reads by their gene id
# -a .gtf: annotation file (from the reference genome)
# -o: output, text and summary file 
apptainer exec --bind /data/ /containers/apptainer/subread_2.0.1--hed695b0_0.sif featureCounts -T 4 -p -s 2 -t exon -g gene_id -a Mus_musculus.GRCm39.113.gtf -o ${XX}_counts_s2_paired.txt ${XX}_sorted.bam