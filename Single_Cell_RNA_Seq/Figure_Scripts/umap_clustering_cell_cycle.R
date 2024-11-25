library(Seurat)
library(SeuratObject)

# Input: filtered seurat object with cell cycle information
DimPlot(cc_phase_seurat, reduction = "umap", group.by = "Phase",
        pt.size = 0.5,
        cols = c("#201547", "#D8AE47", "#483688"))


# Input: filtered seurat object with cell cycle information
DimPlot(cc_phase_seurat, reduction = "umap", group.by = "Phase",
        split.by = "Sample",
        pt.size = 0.5,
        cols = c("#201547", "#D8AE47", "#483688"))

