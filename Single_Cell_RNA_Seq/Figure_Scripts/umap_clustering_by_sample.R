# Load required libraries
library(Seurat)
library(SeuratObject)

# Run UMAP on the filtered and scaled seurat object
seurat_obj_filtered_scaled <- RunUMAP(seurat_obj_filtered_scaled, dims = 1:10)

# Visualize the clusters
DimPlot(seurat_obj_filtered_scaled, reduction = "umap", group.by = "sample",
                       pt.size = 0.5,
                       cols = c("#0d1b53", "#fe828c"))
