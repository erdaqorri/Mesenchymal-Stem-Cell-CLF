# Title: RNA-Seq analysis of canine Mesenchymal Stem Cells
# Date: 25-04-2023

# Load libraries ----
library(rhdf5)
library(tidyverse)
library(tximport)
library(biomaRt)
library(readr)
library(dplyr)
library(matrixStats)
library(DESeq2)
library(edgeR)
library(limma)
library(plotly)

# Read in study design ----
canine_study_design <- read_table("~/path/to/dir/study_design.txt")
canine_study_design <- as_tibble(canine_study_design)

# Load kallisto abundances files
path <- file.path(canine_study_design$Samples, "abundance.tsv")

# Extract gene-level information from biomaRt ----
canis_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                      dataset = "clfamiliaris_gene_ensembl")

# Get attributes
canis_tx_download <- getBM(attributes = c("ensembl_transcript_id_version",
                                          "ensembl_gene_id",
                                          "external_gene_name"),
                           mart = canis_mart)

# Load the kallisto abundance.tsv files ----
Tx_gene_id <- tximport(path,
                      type = "kallisto",
                      tx2gene = canis_tx_download,
                      txOut = FALSE, # dictates if to continue with transcript level or gene level data
                      countsFromAbundance = "lengthScaledTPM")


# Extract differentially expressed genes list
txi_gene_counts <- Tx_gene_id$counts 
deg_list <- DGEList(txi_gene_counts)

# Convert to counts per million
cpm <- cpm(deg_list)
log2_cpm <- cpm(deg_list, log=TRUE)

# Filtering
group <- as.factor(canine_study_design$Group) 
deg_list$samples$group <- group
retain <- filterByExpr(deg_list, group=group)
dgelist_filtered <- deg_list[retain,, keep.lib.sizes=FALSE]

# Normalizing TMM ----
deg_list_filtered_normalized <- calcNormFactors(deg_list_filtered, method = "TMM")
log2_list_filtered_normalized <- cpm(deg_list_filtered_normalized, log=TRUE)

# Principle Component Analysis (PCA) ----
pca_output <- prcomp(t(log2_list_filtered_normalized), scale. = F, retx = T)
pca_output$x

pc_var <- pca_output$sdev^2
pc_pct <- round(pc_var/sum(pc_var)*100, 1)

pca_output.df <- as_tibble(pca_output$x)

pca_plot <- ggplot(pca_output.df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 6, shape = 16) +
  xlab(paste0("PC1 (", pc_pct[1], "%", ")")) +
  ylab(paste0("PC2 (", pc_pct[2], "%", ")")) +
  labs(title = "PCA plot") +  # coord_fixed() +
  theme_classic() +
  ggtitle("") +
  labs(color = "Group", fill = "Group")

# Design Model matrix ----
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Generate Contrast matrix ----
contrast.matrix <- makeContrasts(P6_P2 = Passage6 - Passage2,
                                 levels = design)

# Modeling the mean-variance trend with voom ----
variance_degelist_filtered_normalized <- voom(deg_list_filtered_normalized, design, plot = TRUE)

# Fit a linear model to the data and run differential gene expression analysis ----
fit_linear_model <- lmFit(variance_degelist_filtered_normalized, design)
fits <- contrasts.fit(fit_linear_model, contrast.matrix)
ebfit <- eBayes(fits)
# plotSA(ebfit)

# Extract up-regulated and down-regulated genes
deg_results <- topTable(ebfit,
                        n=Inf)

degs <- subset(deg_results, adj.P.Val < 0.01 & 
                 abs(logFC) > 1)

## Extract up-regulated genes
up_genes <- degs %>% 
  dplyr::filter(logFC > 0)

## Extract down-regulated genes
down_genes <- degs %>%
  dplyr::filter(logFC < 0)