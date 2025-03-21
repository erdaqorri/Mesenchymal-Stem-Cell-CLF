# Code to reproduce the results of our analysis

---

This repository contains the code for the project *Comprehensive Bulk and Single-Cell RNA Sequencing Uncovers Senescence-Associated Biomarkers in Therapeutic Mesenchymal Stem Cells.* 

It includes scripts for:  
- Processing RNA and Single-Cell RNA sequencing data.  
- Data preprocessing and visualization.

## Genome version
Both analyses were conducted using the **Canis lupus familiaris Cfam3.1** transcriptome/genome (updated on November 19, 2019), retrieved from the Ensembl database.

## Single-Cell RNA-Seq Analysis 
The main analysis steps are compiled in an R Markown (`.Rmd`) file.

### Prerequisites
- **CellRanger Version**:  
  The raw FASTQ files were processed using **CellRanger version v.7.2.0**.

- **Seurat Version**:  
  The R Markdown code strictly requires **Seurat version 4.4.0** to compile successfully.

#### Quick Access for Seurat Installation

To install Seurat and its dependencies, you can use the following commands:

```r
remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
```
### Input Data
To reproduce our results, the filtered matrix files are required as input. Replace the placeholder directory paths in the code with the path to your local directories.

### Directory Structure

### Single cell RNA analysis/:
- **`single_cell_analysis_main.Rmd`**: The main analysis R Markdown file contains the steps to filter and process the single-cell RNA-seq data.
- **`cell_cycle_analysis_seurat.R`**: Calculates cell cycle scores used in the cell cycle analysis.
- **`clf_extract_cell_cycle_genes.R`**: Extracts canine-specific cell-cycle genes.

## RNA-Seq Analysis Directory

This directory contains the code to reproduce the RNA-Seq data analysis. All the necessary packages used are included in the R scripts. To run the script, simply replace the paths to your local directories in the code as needed.

### Directory Structure
- **`Preprocessing/`**: Contains bash scripts for quality control and pseudo-alignment.
- **`DGE_Analysis/`**: Contains scripts for differential gene expression analysis, GO/KEGG enrichment, and data visualization.

## Data Availability

1. The raw *Single Cell RNA-Seq* FASTQ files are available in the SRA under accession **PRJNA1235986** and can be downloaded from **[enter link]**.
2. The raw *RNA-Seq* FASTQ files are available in the SRA under accession **PRJNA1235683** and can be downloaded from **[enter link]**.
3. Canis lupus familiaris cell cycle genes are available for download from the Zenodo link: **[enter link]**.

## Contact Information

For any additional information or questions, please feel free to contact Erda Alexandria Qorri at **qorri.erda@brc.hu**.


