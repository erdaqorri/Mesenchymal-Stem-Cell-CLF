# Figure X: Heatmap of the expression of different subclusters
# Date: 26-11-2024

# Define cluster order based on their P2 and P6 origin
cluster_order <- c("0", "3", "4", "7", "1", "2", "5", "6")

# Genes identified using FindAllMarkers function from Seurat, only top 5 were retained
genes <- c("GPC6", "PTX3", "PTGES", "CD44", "GFPT2", "IGFBP2", "TIMP1", "CCND1", "NDUFA4L2", "HSPB1", "COL11A1", "EDIL3", "DIAPH2", "PDLIM5", "PAPPA", "CDC20", "DIAPH3", "STMN1", "RRM2", "HMGB2", "KRT3", "TGM2", "PLA2G7", "IL1RL1", "SFRP2", "DPYSL2", "SEC61G", "FN1", "CEMIP", "PTPN22", "KRT3", "IL1RL1", "S100A5", "CDK14", "OGN", "IGFBP7", "FRZB", "INPP4B", "CRYAB", "TINAGL1")
genes <- unique(toupper(genes)) 

# Extract the expression matrix for the genes of interest
expression_data <- FetchData(filtered_seurat_99_adjusted_gf, , vars = genes)

# Assign cluster information to the extracted expression data
expression_data$Cluster <- filtered_seurat_99_adjusted_gf@meta.data$seurat_clusters

# Aggregate the gene information per cluster and per gene using median
agg_data <- aggregate(. ~ Cluster, data = expression_data, FUN = median)

rownames(agg_data) <- agg_data$Cluster
agg_data$Cluster <- NULL

# Convert intro a matrix for the pheatmap
agg_data_matrix <- as.matrix(agg_data)

# Reorder the matrix rows according to desired cluster order
agg_data_matrix <- agg_data_matrix[cluster_order,]

# Create a row annotation dataframe with ordered factors
row_annotation <- data.frame(
  Cluster = factor(rownames(agg_data_matrix), levels = cluster_order)
)
rownames(row_annotation) <- rownames(agg_data_matrix)

# Define colors for clusters
cluster_colors <- colorRampPalette(c("#4040F0", "#E040E0"))(length(cluster_order))
names(cluster_colors) <- cluster_order

# Create annotation colors list
ann_colors <- list(
  Cluster = cluster_colors
)

# Create the heatmap
pheatmap(t(scale(t(agg_data_matrix))),  # Scale the data
         color = colorRampPalette(c("white", "navy", "red"))(100),
         cluster_rows = FALSE,  # Disable row clustering to maintain our custom order
         cluster_cols = FALSE,  # Don't cluster columns to maintain gene order
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_row = row_annotation,
         annotation_colors = ann_colors,
         fontsize = 10,
         angle_col = 45,
         fontsize_row = 10,
         fontsize_col = 8,
         cellwidth = 15,
         cellheight = 12,
         filename = "heatmap_clusters_median.png",
         width = 12,
         height = 8,
         res = 1200)
