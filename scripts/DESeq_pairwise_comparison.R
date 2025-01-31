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
library(ggplot2)
#------------- Importing and cleaning the count table --------------------------

# Import the count table outputed by featureCounts
count_raw <- read.table("BAM_merged_counts_s2_paired_Q10.txt",
                        header = TRUE, 
                        row.names = 1)

# Trim the first 5 columns contianing information on the genomic features
counts <- count_raw[, -c(1:5)]

#----------- Importing the metadata and renaming columns -----------------------

# The metadata file must contain the sample information in 
# the same order than the count table
metadata <- read.table("metadata.txt", 
                       header = TRUE, 
                       row.names = 1)

# Before the next step, ensure that the metadata contains the sample names 
# in the proper order
colnames(counts) <- rownames(metadata)

# Changing type of experimental conditions to factors for downstream analysis
metadata$Treatment <- as.factor(metadata$Treatment)

#-------- Creating a DESeq object from the counts table and the metadata -------

dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = metadata, 
                              design = ~ Condition)
dds <- DESeq(dds)
# The following functions and plots will be based on this dds object

#------------------------------ PCA plot ---------------------------------------

# Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Extract PCA data and variance explained using plotPCA
pca_data <- plotPCA(vsd, intgroup = c("Tissue", "Treatment"), returnData = TRUE)

percent_var <- round(100 * attr(pca_data, "percentVar"))

# Create the PCA plot with variance in axis labels
ggplot(pca_data, aes(x = PC1, y = PC2, color = Tissue, shape = Treatment)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of infected and control blood and lung samples",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  coord_fixed() +
  theme_minimal()

#------------------------- Defining functions ----------------------------------

# Below are several functions used for plots and DEG analyses

# To perform the pairwise analysis, print the number of DEGs 
# and store the results in a data frame:
pairwise_analysis <- function(contrast, comparison_name, padj_threshold=0.05){
  
  # Comparison based on the groups in contrast
  res_comparison <- results(dds, 
                 contrast = contrast)
  # Generate MA plot
  plotMA(res_comparison, 
         main = paste("MA Plot -", comparison_name))
  
  # Putting the results in a data frame and adding the gene names
  result_dataframe <- add_gene_names(res_comparison)
  
  # Calculate the number of DEGs
  print(comparison_name)
  num_genes <- sum(res_comparison$padj < padj_threshold, na.rm = TRUE)
  print(paste("Total number of differentially expressed genes:", num_genes))
  print(summary(res_comparison, alpha = padj_threshold))
  
  return(result_dataframe)
}

# To retrieve the gene name from Ensembl IDs
add_gene_names <- function(res_comparison) {
  gene_names <- mapIds(
    org.Mm.eg.db,
    keys = rownames(res_comparison),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  # Adding the column GeneName
  res_comparison <- as.data.frame(res_comparison)
  res_comparison$GeneName <- gene_names
  return(res_comparison)
}

# To produce a volcano plot
generating_volcano_plot <- function(dataframe, comparison_name, padj_threshold = 0.05, logfc_threshold = 1) {
  
  # Filter out rows with NAs in log2FoldChange or in adjusted p-values columns
  dataframe_clean <- dataframe %>%
    filter(!is.na(log2FoldChange) & !is.na(padj))
  
  # Create a new column for classification of DEGs (Upregulated, Downregulated, NS)
  dataframe_clean$significance <- ifelse(
    dataframe_clean$padj < padj_threshold & abs(dataframe_clean$log2FoldChange) > logfc_threshold,
    ifelse(dataframe_clean$log2FoldChange > 0, "Upregulated", "Downregulated"), "NS")
  
  # Add a new column for -log10(padj)
  dataframe_clean$log10Padj <- -log10(dataframe_clean$padj)
  
  # Generate the volcano plot
  volcano_plot <- ggplot(dataframe_clean, aes(x = log2FoldChange, y = log10Padj, color = significance)) +
    geom_point(alpha = 0.6, size = 2.5) +
    scale_color_manual(values = c("NS" = "gray", "Upregulated" = "red", "Downregulated" = "blue")) +
    labs(title = paste("Volcano Plot -", comparison_name),
         subtitle = paste("p-value <", padj_threshold, ", |log2FC| >", logfc_threshold),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-Value",
         caption = "Source: DESeq2 Results") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    geom_hline(yintercept = -log10(padj_threshold), 
               linetype = "dashed", 
               color = "black", 
               linewidth = 0.5) +  # p-value threshold line
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dashed", 
               color = "black", 
               linewidth = 0.5) +  # log2FC threshold lines
    guides(color = guide_legend(title = "Gene Regulation"))

  print(volcano_plot)
}

# To print top 10 up and down DEGs: 
top_10_genes <- function(dataframe, comparison_name, logfc_threshold=log(2)) {
  # Filter genes with significant fold change
  upregulated_genes <- dataframe[dataframe$log2FoldChange >= logfc_threshold & dataframe$padj <= 0.05, ]
  downregulated_genes <- dataframe[dataframe$log2FoldChange <= -logfc_threshold & dataframe$padj <= 0.05, ]
  
  # Filter 10 most upregulated genes
  top_upregulated <- upregulated_genes[order(upregulated_genes$log2FoldChange, decreasing = TRUE), ]
  top_upregulated <- head(top_upregulated, 10)
  
  # Filter 10 most downregulated genes 
  top_downregulated <- downregulated_genes[order(downregulated_genes$log2FoldChange), ]
  top_downregulated <- head(top_downregulated, 10)
  
  # Print the results
  print(paste("Top 10 upregulated genes -", comparison_name))
  print(top_upregulated[, c("GeneName", "log2FoldChange", "padj")])
  
  print(paste("Top 10 downregulated genes -", comparison_name))
  print(top_downregulated[, c("GeneName", "log2FoldChange", "padj")])
}

# To perform the GO analysis: 
go_analysis_pairwise <- function(dataframe, comparison_name, ont , regulation="all", output_file=NULL){
  
  ontology_list <- list("BP"="Biological process","CC"="Cellular component","MF"="Molecular function")
  
  # Filter based on regulation type: "up", "down", or "all"
  if (regulation == "up") {
    # Upregulated genes
    genes_de <- rownames(dataframe[dataframe$log2FoldChange >= log(2) & dataframe$padj <= 0.05, ])
  } else if (regulation == "down") {
    # Downregulated genes
    genes_de <- rownames(dataframe[dataframe$log2FoldChange <= -log(2) & dataframe$padj <= 0.05, ])
  } else {
    # All genes with padj <= 0.05
    genes_de <- rownames(dataframe[dataframe$padj <= 0.05 & !is.na(dataframe$padj), ])
  }
  
  # List of all genes measured in the analysis
  genes_universe <- rownames(dataframe)

  # Perform the GO analysis using enrichGO
  go_results <- enrichGO(
    gene = genes_de,
    universe = genes_universe,
    OrgDb = org.Mm.eg.db,
    ont = ont, # "BP", "MF", "CC", or "ALL"
    keyType = "ENSEMBL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2 )
  fold_change <- dataframe$log2FoldChange
  names(fold_change) <- rownames(dataframe)
  
  ontology <- ontology_list[ont]
  # Barplot GO terms
  bar_plot <- barplot(go_results, 
                      title = paste("GO Terms", ontology), 
                      font.size = 12)
  print(bar_plot)
  # Dotplot GO terms
  dot_plot <- dotplot(go_results, 
                      title = paste("GO Terms -", ontology))
  print(dot_plot)
  
  # Get the results for the selected GO category (subontology)
  go_results_df <- go_results@result  # All results for the selected subontology
  
  # Ensure it is a data frame
  go_results_df <- as.data.frame(go_results_df)
  
  # Save results to an Excel file if output_file is specified
  if (!is.null(output_file)) {
    go_results_df <- as.data.frame(go_results@result)
    write_xlsx(list("GO_Results" = go_results_df), path = output_file)
    print(paste("GO results saved to", output_file))
  }
}

# To plot count data for a specific gene as bar plots:
plot_gene_expression <- function(dds_object, gene_name, condition_column = "Condition", keep_levels = NULL) {
  # Ensure the gene exists in the DESeq2 object
  if (!(gene_name %in% rownames(dds_object))) {
    stop(paste("Gene", gene_name, "not found in the DESeq2 object"))
  }
  
  # Create the count plot data for the given gene
  plot_data <- plotCounts(dds_object, gene = gene_name, intgroup = condition_column, returnData = TRUE)
  
  # Filter levels if `keep_levels` is specified
  if (!is.null(keep_levels)) {
    plot_data <- plot_data[plot_data[[condition_column]] %in% keep_levels, ]
    plot_data[[condition_column]] <- droplevels(plot_data[[condition_column]])
  }
  
  # Calculate summary statistics
  summary_data <- plot_data %>%
    group_by(!!sym(condition_column)) %>%
    summarize(mean_count = mean(count), 
              sd_count = sd(count), 
              .groups = "drop")
  
  # Create the plot
  plot <- ggplot() +
    # Add bars for mean counts
    geom_bar(data = summary_data, aes(x = .data[[condition_column]], y = mean_count), 
             stat = "identity", fill = "lightgrey", alpha = 0.8, width = 0.5, color = "black") +
    # Add error bars for standard deviation
    geom_errorbar(data = summary_data, aes(x = .data[[condition_column]], 
                                           ymin = mean_count - sd_count, ymax = mean_count + sd_count),
                  width = 0.2, linewidth = 0.8, color = "black") +
    # Add sample points inside the bars
    geom_point(data = plot_data, aes(x = .data[[condition_column]], y = count), 
               position = position_jitter(width = 0.15, height = 0), size = 3, shape = 21, 
               fill = "white", color = "black", stroke = 0.8) +
    # Log scale for y-axis
    scale_y_log10() +
    # Titles and labels
    ggtitle(paste("Gene Expression for", gene_name)) +
    labs(x = "Experimental Condition", y = "Normalized Counts") +
    # Clean theme
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5), # Straight x-axis labels
      panel.grid.major.x = element_blank(), # Remove vertical grid lines
      panel.grid.minor = element_blank()
    )
  
  # Print the plot
  print(plot)
}

# Convert the Ensemble IDs to gene symbols
prepare_dds_with_symbols <- function(dds) {
  library(org.Mm.eg.db)
  
  # Map ENSEMBL IDs to gene symbols
  gene_symbols <- mapIds(
    org.Mm.eg.db,
    keys = rownames(dds),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Replace NA values with original ENSEMBL IDs
  gene_symbols[is.na(gene_symbols)] <- rownames(dds)[is.na(gene_symbols)]
  
  # Ensure unique gene symbols
  if (any(duplicated(gene_symbols))) {
    gene_symbols[duplicated(gene_symbols)] <- paste0(
      gene_symbols[duplicated(gene_symbols)], "_", rownames(dds)[duplicated(gene_symbols)]
    )
  }
  
  # Assign new rownames to dds
  rownames(dds) <- gene_symbols
  
  # Validate the new dds object
  if (any(is.na(rownames(dds)))) {
    stop("NA values found in rownames of `dds`.")
  }
  if (any(duplicated(rownames(dds)))) {
    stop("Duplicated values found in rownames of `dds`.")
  }
  
  return(dds)
}

#-------------- Pairwise analysis - Lung control vs Lung infected --------------

comparison_name <- "infected lung samples vs control lung sample"

#Produce a specific dataframe for the comparison:
dataframe_lung_controlvsinfected <- pairwise_analysis(c("Condition", "Lung_Infected","Lung_Control"),comparison_name)

#Print top 10 genes up/downregulated
top_10_genes(dataframe_lung_controlvsinfected,comparison_name)

#Produce volcano plot
generating_volcano_plot(dataframe_lung_controlvsinfected,comparison_name)

#Plotting specific gene counts using gene symbols:
dds_gene_symbol <- prepare_dds_with_symbols(dds)
immune_response_genes <- c("Il12a", "Il1a","Il12b", "Il2", "Il2rb", "Ifng", "Tbx21", "Myod1", "Myog", "Gbp3", "Gbp7", "Iigp1", "Tgtp1", "Stat4", "Myh11", "Des", "Pax7")
for (gene in immune_response_genes) {
  plot_gene_expression(dds_gene_symbol, gene, keep_levels = c("Lung_Control", "Lung_Infected"))
}
results_object <- results(dds_gene_symbol, contrast = c("Condition", "Lung_Infected", "Lung_Control"))
specific_genes <- results_object[rownames(results_object) %in% immune_response_genes, ]
specific_genes

#Go analysis: barplots, dotplots and excel reports
go_analysis_pairwise(dataframe_lung_controlvsinfected,comparison_name,"BP", output_file = "GO_Enrichment_Lung_BP.xlsx")
go_analysis_pairwise(dataframe_lung_controlvsinfected,comparison_name,"MF", output_file = "GO_Enrichment_Lung_MF.xlsx")
go_analysis_pairwise(dataframe_lung_controlvsinfected,comparison_name,"CC", output_file = "GO_Enrichment_Lung_CC.xlsx")

#------------- Pairwise analysis - Blood control vs Blood infected -------------

comparison_name <- "infected blood samples vs control blood sample"

#Produce a specific dataframe for the comparison:
dataframe_blood_controlvsinfected <- pairwise_analysis(c("Condition", "Blood_Infected","Blood_Control"),comparison_name)
#Print top 10 genes up/down-regulated
top_10_genes(dataframe_blood_controlvsinfected,comparison_name)
#Produce volcano plot
generating_volcano_plot(dataframe_blood_controlvsinfected,comparison_name)

#Go analysis
go_analysis_pairwise(dataframe_blood_controlvsinfected,comparison_name,"BP", output_file = "GO_Enrichment_Blood_BP.xlsx")
go_analysis_pairwise(dataframe_blood_controlvsinfected,comparison_name,"MF", output_file = "GO_Enrichment_Blood_MF.xlsx") 
go_analysis_pairwise(dataframe_blood_controlvsinfected,comparison_name,"CC", output_file = "GO_Enrichment_Blood_CC.xlsx") 

#--------------- Pairwise analysis - Blood control vs Lung control -------------
comparaison_name = "blood vs lung control"
#Produce a specific dataframe for the comparison:
dataframe_control_bloodvslung <- run_pairwise_analysis(c("Condition", "Blood_Control","Lung_Control"),comparison_name)
generating_volcano_plot(dataframe_control_bloodvslung,comparison_name)

#--------------- Pairwise analysis - Blood infected vs Lung infected -----------
comparaison_name = "blood vs lung infected"
#Produce a specific dataframe for the comparison:
dataframe_infected_bloodvslung <- run_pairwise_analysis(c("Condition", "Blood_Infected","Lung_Infected"),comparison_name)
generating_volcano_plot(dataframe_control_bloodvslung,comparison_name)
