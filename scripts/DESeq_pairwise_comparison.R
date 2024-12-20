#-------------------------------------------------------------------------------
#-------------- DESeq2 Differential Expression Analysis ------------------------
#-------------------------------------------------------------------------------

#---------------------- Loading libraries --------------------------------------
library(clusterProfiler)
library(pheatmap)
library(stringr)
library(DESeq2)
library(RColorBrewer)
library(EnhancedVolcano)
library(org.Mm.eg.db)
library(tidyverse)
library(writexl)

#----------- Importing and cleaning the count table ----------------------------

count_raw <- read.table("BAM_merged_counts_s2_paired_Q10.txt",
                        header = TRUE, 
                        row.names = 1)

#There are 5 columns at the beginning that are not of interest for us
counts <- count_raw[, -c(1:5)]

#--------- Importing the metadata and renaming columns -------------------------

metadata <- read.table("metadata.txt", 
                       header = TRUE, 
                       row.names = 1)

#The metadata file must contain the sample information in the same order than the count table
#Before this step, ensure that the metadata contains the sample names in the proper order
colnames(counts) <- rownames(metadata)

#Changing type of experimental conditions to factors for downstream analysis
metadata$Condition <- as.factor(metadata$Condition)
metadata$Tissue <- as.factor(metadata$Tissue)

#-------- Creating a DESeq object from the counts table and the metadata -------

dds_Cond_Tissue_All <- DESeqDataSetFromMatrix(countData = counts, 
                                              colData = metadata, 
                                              design = ~ Tissue * Condition )
#Performing the DESeq analysis
dds_Cond_Tissue_All <- DESeq(dds_Cond_Tissue_All)

#------ Pairwise comparison of infected vs control lung samples ----------------

#Select only the data from Lung samples
dds_lung <- dds_Cond_Tissue_All[, metadata$Tissue == "Lung"]

#Proceed to another DESeq analysis, but on this subsample only
dds_lung <- DESeqDataSetFromMatrix(
  countData = counts(dds_lung),
  colData = colData(dds_lung),
  design = ~ Condition)

dds_lung <- DESeq(dds_lung)

#Compute the result table for the pairwise comparison between lung samples
res_lung <- results(dds_lung, contrast = c("Condition", "Infected", "Control"))

#Number of significant genes at padj < 0.05:
sum(res_lung$padj < 0.05, na.rm = TRUE)

#Log the number of genes that are up or downregualted
upregulated <- sum(res_lung$padj < 0.05 & res_lung$log2FoldChange > 0, na.rm = TRUE)
downregulated <- sum(res_lung$padj < 0.05 & res_lung$log2FoldChange < 0, na.rm = TRUE)

cat("Upregulated:", upregulated, "\nDownregulated:", downregulated)

#-------------- Changing Ensembl IDs with gene symbols -------------------------

#Extract Ensembl IDs from the rownames of the DESeq object
ens <- rownames(dds_lung)

#Map Ensembl IDs to gene symbols
symbols <- mapIds(org.Mm.eg.db, 
                  keys = ens, 
                  column = "SYMBOL", 
                  keytype = "ENSEMBL")

#Log the number of unmapped Ensembl IDs
cat("Number of unmapped Ensembl IDs:", sum(is.na(symbols)), "\n")

#Replace NA values with the original Ensembl IDs for genes that couldn't be mapped
symbols[is.na(symbols)] <- rownames(dds_lung)[is.na(symbols)]

#Reorder symbols to match the original data
symbols <- symbols[match(rownames(dds_lung), names(symbols))]

#Replace rownames of dds_lung with gene symbols
rownames(dds_lung) <- symbols

#Verify the updated rownames
head(rownames(dds_lung))

#------------------------- Visualize specific genes ----------------------------

#Extract the top 5 upregulated genes (sorted by log2FoldChange)
top_upregulated_genes <- rownames(res_lung[order(res_lung$log2FoldChange, decreasing = TRUE), ])[1:5]

#Extract the top 5 downregulated genes (sorted by log2FoldChange)
top_downregulated_genes <- rownames(res_lung[order(res_lung$log2FoldChange), ])[1:5]

#Combine into one vector
top_genes <- c(top_upregulated_genes, top_downregulated_genes)

#Log the selected top genes
print(top_genes)

#Create count plots for these top genes:
for (gene in top_genes) {
  #Create a count plot for the current gene
  plot <- plotCounts(dds_lung, gene = gene, intgroup = "Condition", returnData = TRUE)
  p <- ggplot(plot, aes(x = Condition, y = count)) + 
    geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) +  #Jitter to avoid overlap
    scale_y_log10(breaks = c(25, 100, 400)) +  #Log scale for counts
    ggtitle(paste("Gene Expression for", gene))
  print(p)  #Print the plot for each gene
}

#-------------------------- Volcano plot ---------------------------------------

#Convert the DESeq2 result to a data frame
volcano_data <- as.data.frame(res_lung)

#Filter out rows with NA in log2FoldChange or padj columns
volcano_data_clean <- volcano_data %>%
  filter(!is.na(log2FoldChange) & !is.na(padj))
head(volcano_data_clean)

#Create a new column for classification of DEGs
volcano_data_clean$significance <- ifelse(
  volcano_data_clean$padj < 0.05 & abs(volcano_data_clean$log2FoldChange) > 1,
  ifelse(volcano_data_clean$log2FoldChange > 0, "Upregulated", "Downregulated"),
  "NS"
)

#Add a new column for -log10(padj)
volcano_data_clean$log10Padj <- -log10(volcano_data_clean$padj)

#Create a volcano plot
ggplot(volcano_data_clean, aes(x = log2FoldChange, y = log10Padj, color = significance)) +
  geom_point(alpha = 0.6, size = 2.5) +  #Adjust point transparency and size
  scale_color_manual(values = c("gray" = "gray", "Upregulated" = "red", "Downregulated" = "blue")) +
  labs(
    title = 'Volcano Plot of DEGs in Infected vs Control Lungs',
    subtitle = 'p-value < 0.05, |log2FC| > 1',
    x = 'Log2 Fold Change',
    y = '-Log10 Adjusted P-Value',
    caption = 'Source: DESeq2 Results'
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +  #p-value threshold line
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.5) +  #log2FC threshold lines
  guides(color = guide_legend(title = "Gene Regulation"))

#---------------------- Functional analysis on selected DEGs -------------------

#Create the list of selected DEGs
res_lung_df <- as.data.frame(res_lung)
# Filter rows where padj < 0.05 and padj is not NA
significant_genes_dds_lung <- res_lung_df[!is.na(res_lung_df$padj) & res_lung_df$padj < 0.05, ]


#Gene identifiers are in the first column of the dataframe
gene_identifiers <- rownames(significant_genes_dds_lung)

#Separate Ensembl IDs and gene symbols
ensembl_ids <- gene_identifiers[grep("^ENSMUSG", gene_identifiers)]  #IDs starting with 'ENSMUSG'
gene_symbols <- gene_identifiers[!gene_identifiers %in% ensembl_ids] #The rest assumed to be symbols

#Log how many IDs were separated
cat("Number of Ensembl IDs:", length(ensembl_ids), "\n")
cat("Number of Gene Symbols:", length(gene_symbols), "\n")

#Map gene symbols to Entrez IDs
entrez_from_symbols <- mapIds(org.Mm.eg.db, 
                              keys = gene_symbols, 
                              column = "ENTREZID", 
                              keytype = "SYMBOL", 
                              multiVals = "first")

cat("Number of genes successfully mapped to Entrez IDs:", sum(!is.na(entrez_from_symbols)), "\n")
#Filter out any NA values (genes that could not be mapped)
entrez_from_symbols_clean <- entrez_from_symbols[!is.na(entrez_from_symbols)]

#Log the number of successfully mapped genes
cat("Number of genes successfully mapped to Entrez IDs:", length(entrez_from_symbols_clean), "\n")


#Perform GO enrichment for BP, MF and CC
go_results_bp <- enrichGO(
  gene          = entrez_from_symbols_clean,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

go_results_mf <- enrichGO(
  gene          = entrez_from_symbols_clean,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

go_results_cc <- enrichGO(
  gene          = entrez_from_symbols_clean,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

#Create individual plots for each subontology
dotplot(go_results_bp, showCategory = 10)
dotplot(go_results_mf, showCategory = 10)
dotplot(go_results_cc, showCategory = 10)

barplot(go_results_bp, showCategory = 10)
barplot(go_results_mf, showCategory = 10)
barplot(go_results_cc, showCategory = 10)

#Finally output the results

#Get the results for each GO category
go_results_bp_df <- go_results_bp@result  #BP (Biological Process)
go_results_mf_df <- go_results_mf@result  #MF (Molecular Function)
go_results_cc_df <- go_results_cc@result  #CC (Cellular Component)

go_results_bp_df <- as.data.frame(go_results_bp_df)
go_results_mf_df <- as.data.frame(go_results_mf_df)
go_results_cc_df <- as.data.frame(go_results_cc_df)


#Create a list of results
results_list <- list(
  "GO_Biological_Process" = go_results_bp_df,
  "GO_Molecular_Function" = go_results_mf_df,
  "GO_Cellular_Component" = go_results_cc_df
)

#Write all results to an Excel file with multiple sheets
write_xlsx(results_list, path = "GO_Enrichment_Results_Lung_CORRECTED.xlsx")



