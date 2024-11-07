# Detect differentially expressed genes from bulk RNA-seq data
Workflow of differential expression analysis on mice infected by Toxoplasma gondii

# Table describing experimental groups for each sample. All samples are from wildtype (WT) mice. 
Since the reads were produced in paired-end mode, there are 2 files per sample, with read 1 and read 2 respectively.

Sample = ID as it is in the fastq file name
Case = Infected
Control = Uninfected


Sample	Group
SRR7821921	Lung_WT_Case
SRR7821922	Lung_WT_Case
SRR7821918	Lung_WT_Case
SRR7821919	Lung_WT_Case
SRR7821920	Lung_WT_Case
SRR7821937	Lung_WT_Control
SRR7821938	Lung_WT_Control
SRR7821939	Lung_WT_Control
SRR7821949	Blood_WT_Case
SRR7821950	Blood_WT_Case
SRR7821951	Blood_WT_Case
SRR7821952	Blood_WT_Case
SRR7821953	Blood_WT_Case
SRR7821968	Blood_WT_Control
SRR7821969	Blood_WT_Control
SRR7821970	Blood_WT_Control


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