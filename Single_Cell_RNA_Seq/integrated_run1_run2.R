####################################### Run 1 and Run 2 Integration ####################################### 

library(SingleCellExperiment)
library(Seurat)
library(SeuratObject)
library(tidyverse) 
library(Matrix)
library(scales)
library(cowplot) 
library(RCurl) 
library(magrittr) 
library(ggplot2)
library(dplyr) 
library(patchwork) 
library(ggpubr)
library(glmGamPoi)

# Load the filtered seurat object assays for both the first and second run
run_1_filtered <- readRDS("filtered_run1.rsd")
run_2_no_mito_filtered <- readRDS("filtered_run2.rds")

# The filtered matrices from the second run do not contain the mitochondrial reads and information which was taken into account in the first run
# The reads mapping to the mitochondrial genes were removed from first run to avoid bias during feature gene and anchor gene identification

# Removing MT genes from the recent dataset ----
dog_mt_genes <- c("MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6",
                  "COX3", "MT-ND3", "ND4L", "MT-ND4", "ND5", "MT-ND6", "MT-CYB")

mt_in_data <- intersect(dog_mt_genes, rownames(run_1_filtered))

run_1_no_mito_filtered <- subset(
  run_1_filtered,
  features = setdiff(rownames(run_1_filtered), mt_in_data)
)

# Check if the filtering worked properly
AverageExpression(run_1_no_mito_filtered, features = dog_mt_genes)$RNA

#################################### Integration with no mt genes ####################################
seurat_obj_list <- list("Run1" = run_1_no_mito_filtered, "Run2" = run_2_no_mito_filtered)

# Assign identity to the cells based on the run
merged_obj <- merge(
  x = run_1_no_mito_filtered,
  y = run_2_no_mito_filtered,
  add.cell.ids = c("Run_1", "Run_2") # optional: prefix cell names
)

merged_obj@meta.data <- merged_obj@meta.data %>%
  mutate(Run_Number = case_when(
    str_detect(sample, "^P372") ~ "P2-R2",
    str_detect(sample, "^P376") ~ "P6-R2",
    str_detect(sample, "^P2") ~ "P2-R1",
    str_detect(sample, "^P6") ~ "P6-R1",
    TRUE ~ NA_character_
  ))

# Split the seurat assays to apply the SCTransform separately
split_seurat <- SplitObject(merged_obj, split.by = "Run_Number")

split_seurat_sct <- lapply(X = split_seurat, FUN = function(x) {
  x <- SCTransform(x, variable.features.n = 3000)
})

# Identify integration features, set to 3000 features
integration_features <- SelectIntegrationFeatures(object.list = split_seurat_sct, nfeatures = 3000)

# Prepare the seurat object for integration
split_seurat_sct <- PrepSCTIntegration(object.list = split_seurat_sct, anchor.features = integration_features)

integration_anchors <- FindIntegrationAnchors(object.list = split_seurat_sct, 
                                         normalization.method = "SCT",
                                         anchor.features = integration_features)

# Integrate the data using the identified anchor markers
integrated_seurat_sct <- IntegrateData(anchorset = integration_anchors, normalization.method = "SCT")
DimPlot(integrated_run,
        split.by = "sample") 

# Set the default assay to integrated
DefaultAssay(integrated_seurat_sct) <- "integrated"

# Determining the best PCs for integration
# Set seed
set.seed(123456)

# Apply the standard pipeline
integrated_seurat_sct <- ScaleData(integrated_seurat_sct, verbose = FALSE)

# PC30 allows for good per sample integration and manages to give high subcluster resolution
integrated_seurat_sct <- RunPCA(integrated_seurat_sct, npcs = 30, verbose = FALSE)
integrated_seurat_sct <- RunUMAP(integrated_seurat_sct, reduction = "pca", dims = 1:30)
integrated_seurat_sct <- FindNeighbors(integrated_seurat_sct, reduction = "pca", dims = 1:30)
integrated_seurat_sct <- FindClusters(integrated_seurat_sct, resolution = 0.2)

# Visualize UMAP
DimPlot(integrated_run, reduction = "umap", group.by = "Run_Number",
                              pt.size = 0.5,
        cols = c("#FFD700", "#5c7650", "green", "#e94560"))

# Visualize the subclusters UMAP
DimPlot(integrated_run, reduction = "umap",
                          label.box = TRUE,
                          label = TRUE,
                          label.size = 4)


# Visualize the number of cells per subcluster to check that they are contributing equally to the separation
ggplot(integrated_run@meta.data) +
  geom_bar(aes(x=integrated_snn_res.0.2, fill=Run_Number), position=position_fill()) +
  theme_classic() +
  labs(x = "Cluster Identity",
       y = "Cell Proportion",
       fill = "Sample ID") +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)
  )

# Check the default assay, for visualization the RNA assay normalized is used
DefaultAssay(integrated_run)

# Reset to RNA if integrated or SCT
DefaultAssay(integrated_run) <- "RNA"

# Normalize for visualization purposes
integrated_run <- NormalizeData(integrated_run, verbose = FALSE)
integrated_run <- ScaleData(integrated_run)

# Visualize the markers
FeaturePlot(integrated_run, features = markers_part1)
FeaturePlot(integrated_run, features = markers_part2)

# Identify markers between the P6 and P2 subclusters
Idents(integrated_run) <- integrated_run$integrated_snn_res.0.2

cluster.P2.P6.markers <- FindMarkers(integrated_run, 
                                  ident.1 = c(0, 4, 5, 6, 7), 
                                  ident.2 = c(1, 2, 3, 8),
                                  only.pos = FALSE,
                                  logfc.threshold = 0.25,
                                  min.pct = 0.25)
