# Load required libraries
library(Seurat)
library(SeuratObject)
library(patchwork)

# Identify the cell populations using resolution 0.2
seurat_obj_filtered_scaled <- FindClusters(seurat_obj_filtered_scaled, resolution = c(0.2))

Idents(object = seurat_obj_filtered_scaled) <- "RNA_snn_res.0.2"

# Plot
DimPlot(seurat_obj_filtered_scaled, reduction = "umap", label = TRUE,
        label.size = 7,
        pt.size = 0.5,
        label.color = "black",
        cols = c("#ad497b", "#cb98dc", "#a6cd8e", "#3496BB", "#92e3de", "#5c7650", "#e94560", "#FFFF00FD"))
