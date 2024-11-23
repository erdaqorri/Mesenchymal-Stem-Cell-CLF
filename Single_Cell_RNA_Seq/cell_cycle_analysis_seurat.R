library(readr)
library(readxl)
library(tidyverse)
library(dplyr)
library(AnnotationDbi)
library(AnnotationHub)
library(org.Cf.eg.db)
library(Seurat)

# Read in the cell cycle genes identified with a custom pipeline
cell_cycle_genes_dog <- read_excel("~/Downloads/clfamiliaris_cell_cycle_genes_25-10-2022.xlsx")

cell_cycle_genes_dog$Gene_Symbol <- mapIds(x = org.Cf.eg.db,
                                           keys = cell_cycle_genes_dog$geneID,
                                           column = "SYMBOL",
                                           keytype = "ENSEMBL",
                                           multiVals = "first")

# Retrieve G2M phase genes
g2m_genes <- cell_cycle_genes_dog %>% dplyr::filter(phase == "G2/M") %>%
  dplyr::select(c(phase, Gene_Symbol, geneID))

# Retrieve S phase genes
s_genes <- cell_cycle_genes_dog %>% dplyr::filter(phase == "S") %>%
  dplyr::select(c(phase, Gene_Symbol, geneID))

if (nrow(g2m_genes) ==0 | nrow(s_genes) == 0) {
  stop("G2M or/and S gene files are empty. Check the input file!")
}

# Assign Cell Cycle scores to cells
cc_phase_seurat <- CellCycleScoring(filtered_seurat_99_adjusted_gf, 
                                    g2m.features = g2m_genes$Gene_Symbol, 
                                    s.features = s_genes$Gene_Symbol,
                                    set.ident = TRUE)
