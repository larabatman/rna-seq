#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=03:00:00
#SBATCH --partition=pibu_el8


#It appears that there is a small typo in the container pathway:
#/containers/apptainer/subread_2.0.1--hed695b0_0.sif

#This version of featureCounts does not have the following option
#featureCounts -p --countReadPairs -t exon -g gene_id -a annotation.gtf -o counts.txt mapping_results_PE.bam


#Directory of the BAM files
BAM_DIR="/data/users/lland/rna_seq/ref_genome/mapped_reads"
cd ${BAM_DIR}

#Create a variable containing the names of all the sorted BAM file to proceed with a merged count
BAM_FILES=$(ls *_sorted.bam)


# Run featureCounts with the following options:
# -T 4: four threads
# -p: paired-end
# -s 2: strandedness
# -t exon: -t to define the feature to count; here, we want to count exons
# -g gene_id: to group the reads by their gene id
# -a .gtf: annotation file (from the reference genome)
# -o: output, text and summary file 
#optionally, we can add a filtering step with -Q 10 that will filter out the reads for which the quality is below 10
apptainer exec --bind /data/ /containers/apptainer/subread_2.0.1--hed695b0_0.sif featureCounts -T 4 -p -s 2 -Q 10 -t exon -g gene_id -a Mus_musculus.GRCm39.113.gtf -o BAM_merged_counts_s2_paired_Q10.txt ${BAM_FILES}
