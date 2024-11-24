# Title: Generate volcano plot
# Date: 25-04-2023

# Load libraries ----
library(EnhancedVolcano)
library(dplyr)
library(tidyverse)
library(limma)
library(magrittr)

# Read in the ebfit object ----
ebfit <- readRDS(file = "/path/to/dir/ebfit_object.rds")

# Sort genes by logFC ----
every_gene <- topTable(ebfit, adjust.method = "BH", coef = 1, number = 20000, sort.by = "logFC")
every_gene_df <- as.data.frame(every_gene)
every_gene_df <- cbind(geneID = rownames(every_gene_df), every_gene_df)
rownames(every_gene_df) <- NULL

# Generate volcano plot ----
EnhancedVolcano(every_gene_df,
                lab = every.gene.df$geneID,
                x = 'logFC',
                y = 'adj.P.Val',
                title = "P6 vs P2",
                subtitle = "",
                pCutoff = 0.01,
                FCcutoff = 1,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendPosition = "right",
                legendLabSize = 12,
                legendIconSize = 4.0)
