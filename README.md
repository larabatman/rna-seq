# Detect differentially expressed genes from bulk RNA-seq data
Workflow of differential expression analysis on mice infected by Toxoplasma gondii.
The workflow is mainly based on the RNA-seq course from the university of Bern (course code 467713), following the steps described in file:///Users/larabatman/Desktop/UNIFR/FS12024/RNA_seq/DEAnalysisWorkflow_toxoplasma.html


## Table of Content

0. [Samples and experimental conditions](#samples-and-experimental-conditions)
1. [Quality Control](#quality-control)
2. [Mapping to the reference genome](#mapping-to-the-reference-genome)
3. [Counting the number of reads per gene](#counting-the-number-of-reads-per-gene)
4. [Differential Expression Analysis](#differential-expression-analysis)

## Samples and experimental conditions
All samples are from wildtype  (WT) mice. 
Since the reads were produced in paired-end mode, there are 2 files per sample, with read 1 and read 2 respectively.


| Sample    | Tissue  | Condition | Replicate |
|-----------|---------|-----------|-----------|
| SRR7821918 | Lung    | Infected  | 1         |
| SRR7821919 | Lung    | Infected  | 2         |
| SRR7821920 | Lung    | Infected  | 3         |
| SRR7821921 | Lung    | Infected  | 4         |
| SRR7821922 | Lung    | Infected  | 5         |
| SRR7821937 | Lung    | Control   | 1         |
| SRR7821938 | Lung    | Control   | 2         |
| SRR7821939 | Lung    | Control   | 3         |
| SRR7821949 | Blood   | Infected  | 1         |
| SRR7821950 | Blood   | Infected  | 2         |
| SRR7821951 | Blood   | Infected  | 3         |
| SRR7821952 | Blood   | Infected  | 4         |
| SRR7821953 | Blood   | Infected  | 5         |
| SRR7821968 | Blood   | Control   | 1         |
| SRR7821969 | Blood   | Control   | 2         |
| SRR7821970 | Blood   | Control   | 3         |

Table 1: Experimental conditions and sample identifiers.

## Quality Control

Let's start by assessing the quality of our raw reads with FastQC.

Run script 01a_QC_raw_reads.sh

Once the fastaqc report is produced for all the samples, it is possible to run a MultiQC analysis in order to visualize more easily the overall quality of our reads.

Run script 01b_MultiQC_raw_reads.sh

Upon inspection, it appears that the overall quality of the reads is good but there are a lot of duplicated reads. As the per base sequenc content suggest, most sequences have discrepancies in their bases percentage at the beginning and end of the reads. It could be possible theoretically to trim the beginning and end of the reads. That has been explored by running file 01c_trimming_raw_reads.sh and 01d_MultiQC_trimmed_reads.sh, using Trimmomatic with different parameters, but no significant difference was noted in the MultiQC report between trimmed and non-trimmed reads. 
The duplicated reads can be managed after alignement to the reference genome, to assess their truly duplicated nature. 

## Mapping to the reference genome

After quality control ensured the overall good quality of our sequenced data, we can align the reads of each sample to the reference genome of Mus musculus. The strain of mice used for sequencing was C57BL/6J. However, we are using the reference genome from the Ensembl database, which contains the assembly DNA for GRCm39. Thus, the overall alignement rate might encount for those differences. 
The first step is to download the latest available versio of the reference genome, both fasta and gtf files. That can be achived through the following commands: 
wget https://ftp.ensembl.org/pub/release-113/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz

Alternatively, it is possible to run the 00_Reference_Genome.sh file that will retrieve and unzip the fasta and GTF files for Mus musculus form the Ensembl database. 

HiSAT2 tools will be used for indexing of he reference genome as well as alignement of the sequenced samples to the reference genome.

To build the reference genome index: 
Run script 02a_Indexing_ref_genome.sh

To align the sequenced samples to the reference genome:
Run script 02b_Mapping.sh

The overall alignement rates of the samples to the reference genome goes from 87.90% to 98.09%. 

HiSAT produces SAM files that need to be converted to BAM files. SAMTOOLS view will handle the conversion:

Run 02c_Convert_SAM_to_BAM.sh

Then SAMTOOLS sort is used to arrange the alignements in each bam file according to their genomic position:

Run 02d_Sorting_BAM.sh

With the sorted bam files, we can then create indexes for quick reference to specific regions of the sorted bam files using SAMTOOLS index:

Run 02e_Indexing_BAM.sh

We have now produced an indexed alignement of each sample to the reference genome of Mus musculus. The alignements can be visualized using IGV by dowloading the respective bam files. That would require to download in your local computer the BAM file containing all the mapped reads, which can take a lot of space - we can instead cut a chromosome of one of the samples to check the alignement visually. An index for that chromosome is also required. 

Run 02f_optional_IGV.sh

## Counting the number of reads per gene

Once the BAM files have been sorted and indexed, we can finally count the numbers of reads that aligned to exons in  order to compare counts from different genes and effectively perform differential expression analysis. 
The featureCounts tool will be used to count reads that mapped to exons. This tool needs the GTF file (unzipped) from the reference genome as well as the sorted BAM file for each sequence. As our data is paired-end and reverse stranded, we are using options -p and -s 2 to respectively account for that. FeatureCounts will produce two text files: one containing the counts for each gene, and a summary file specifying the amount of reads that were succesfully aligned and other metrics. 

Run 03_Counts.sh

Actually, featureCounts also offers the possibility to merge samples into one count table, which would be helpful here. 

Run 03a_Counts_all_BAM.sh

Now that the heavy computing part of the workflow is done (producing the count table from the fastq files of each sample), we can move on to the actual differential expression analysis which will be performed using the DESeq2 package in R.

## Differential Expression Analysis

After producing the count table using featureCount, we can finally move on from the cluster to our local computer and perform the differential expression analysis using the DESeq2 R package.
For this, we will need to create a metadata text file that will contain the description of our samples (sample IDs with theri corresponding condition (Infected or Control) and tissue type (Lung or Blood), as well as their replicate number), much like Table 1.
We will also need to copy the counts table from the cluster to our local machine. 
Once we have these two files, we can run the analysis:

On R, run the DESeq_pairwise_comparison.r

The script first imports and trims the count table to ensure that it is in the format accepted by DESeq2. Then, the metadata file is also exported. 
The dds object is created with the DESeqDataSetFromMatrix function, which creates a data set with the count data, the metadata and a design formula that accounts for the experimental conditions. Here, we use the parameter "Condition" which descibres the state (Infected or Control) as well as the tissue (Lung or Blood).
The DESeq function is called on the dds object, such that the differential expression analysis performed: the log2 Fold Changes are generated and their respective adjusted p-values are calculated for each gene across all samples. 
The dds object is also used, after stabilizing its variance, for PCA analysis.
For a given pairwise comparison, for example infected lung samples against control lung samples, a dataframe containing relevant information form the analyzed dds object is produced through the function result and the specified contrast parameters (Lung_Infected and Lung_Condition). The number of DEGs is calculated. The top 10 most upregulated and downregulated genes are then determined, and a volcano plot is produced. Genes of interest can be manually searched, and count barplots can be generated as well. 
Finally, the over-representation of GO terms analysis produces barplots, fot plots and Excel reports for the specified comparison. 

<!--

#Check the status of your repository
git status

#Stage changes
git add filename   # or use 'git add .' to stage all changes

#Commit changes
git commit -m "Your commit message"

#Push to GitHub
git push origin main

#To save changes before pulling
git stash


#To retrieve the changes after pulling
git stash pop
-->

