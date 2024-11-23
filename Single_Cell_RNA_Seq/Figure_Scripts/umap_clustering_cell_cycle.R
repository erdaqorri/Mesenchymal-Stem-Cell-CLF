library(Seurat)
library(SeuratObject)

# Load seurat object stored as RSD
seurat_obj <-  read_rds("/path/to/rsd_file/seurat_obj_cell_cycle.rsd")

# Merged 
DimPlot(cc_phase_seurat, reduction = "umap", group.by = "Phase",
        pt.size = 0.5,
        cols = c("#201547", "#D8AE47", "#483688"))


# Split by samples
DimPlot(cc_phase_seurat, reduction = "umap", group.by = "Phase",
        split.by = "Sample",
        pt.size = 0.5,
        cols = c("#201547", "#D8AE47", "#483688"))

