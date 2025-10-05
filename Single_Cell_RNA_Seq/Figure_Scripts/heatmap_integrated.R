# Check the default assay
DefaultAssay(integrated_seurat_sct)

# Scale the data
integrated_seurat_sct <- ScaleData(integrated_seurat_sct)

# List of upregulated marker genes
upregulated <- c(
  "IGFBP2", "COL11A1", "TIMP1", "CCND1", "CRYAB", "NDUFA4L2", "THBS1", 
  "INHBA", "CDKN1A", "FN1", "EDIL3", "VCAN", "CDKN2A", "TGFB2", "COL8A1", 
  "TPM2", "CCN2", "EDN1", "S100A8", "HSPB1", "HMGB1", "CYP7B1", "SMG6", 
  "PARD3B", "GHR", "RRM2", "DCLK1", "NFATC2", "UBE2S", "SPC24", "MYBL2", 
  "DIAPH3", "BNC2", "CENPW", "STMN1", "PLA2G7", "DHDDS", "GPC6", "PTGES", "CDC20"
)

# Set cluster identities
Idents(integrated_seurat_sct) <- "seurat_clusters"

# Define cluster order (from P2 to P6)
cluster_order <- c(1, 2, 3, 8, 0, 4, 5, 6, 7)

# Apply custom cluster ordering
integrated_seurat_sct$seurat_clusters <- factor(
  integrated_seurat_sct$seurat_clusters,
  levels = cluster_order
)

# Generate heatmap
DoHeatmap(
  integrated_seurat_sct,
  features = upregulated,
  group.by = "seurat_clusters",
  size = 5,
  angle = 0,
  draw.lines = TRUE,
  lines.width = 2,
  group.bar.height = 0.02
) +
  theme(
    axis.text.y = element_text(size = 12, angle = 0, hjust = 1, color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
