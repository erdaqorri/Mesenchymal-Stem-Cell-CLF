# Title: genes2clfamiliaris
# Based on the script from Victor Barrera

# Loading libraries
library(biomaRt)
library(tidyverse)
library(RCurl)
library(dplyr)
library(rio)

## Download latest version of the human cell cycle genes
human_cell_cycle_genes <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv")
human_cell_cycle_genes <- read.csv(text = human_cell_cycle_genes)

## Download Ensembl data for hsapiens
ensembl <- useMart("ensembl") # creates an object of the mart class
ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # probably extracts all the human genes

## Download Ensembl data for clfamiliaris
ensembl.bork <- useMart("ensembl")
ensembl.bork.doggo <- useMart("ensembl", dataset = "clfamiliaris_gene_ensembl")

# Obtaining orthologue
doggo_ortholog_information <-
  getBM(
    attributes = c(
      'ensembl_gene_id',
      'clfamiliaris_homolog_ensembl_gene',
      'clfamiliaris_homolog_orthology_type',
      'clfamiliaris_homolog_orthology_confidence'
    ),
    filters = 'ensembl_gene_id',
    values = human_cell_cycle_genes$geneID,
    mart = ensembl.human
  )

doggo_ortholog_information %>%
  View()

# We only keep the ones with one2one orthology type and high orthology confidence
doggo_selected_genes <-
  doggo_ortholog_information %>% dplyr::filter(clfamiliaris_homolog_orthology_type == "ortholog_one2one",
                                  clfamiliaris_homolog_orthology_confidence == 1)

doggo_selected_genes %>%
  View()

## Adding Cell Cycle information
doggo_cc_genes <- doggo_selected_genes %>% dplyr::inner_join(human_cc_genes, by = c("ensembl_gene_id" = "geneID")) %>% 
  dplyr::select(phase,geneID=clfamiliaris_homolog_ensembl_gene) %>% dplyr::arrange(phase)

doggo_cc_genes$modified <- Sys.Date()

# Save file
rio::export(doggo_cc_genes, file = file.path("/home/user/10X_Genomics_Data/", "clfamiliaris.csv"))
