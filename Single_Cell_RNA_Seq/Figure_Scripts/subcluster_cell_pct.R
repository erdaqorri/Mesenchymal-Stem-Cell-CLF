# Supplementary Figure 4: Distribution of cells across every cell subcluster 

library(Seurat)
library(dplyr)
library(tidyverse)
library(magrittr)
library(ggplot2)

# Read in the filtered and processed seurat object
filtered_seurat_99_adjusted_gf <- read_rds("/path/filtered_seurat_object_99.rsd")

# Extract cluster information
cluster_information <- Idents(filtered_seurat_99_adjusted_gf)

# Extract cells per cluster and sample
samples <- c("P2", "P6")

subcluster_cells_list <- list()

for (sample_id in samples) {
   
   seurat_subset <- subset(filtered_seurat_99_adjusted_gf, subset = sample == sample_id)
   seurat_subcluster_id <- Idents(seurat_subset)
   subcluster_cells_list[[sample_id]] <- as.data.frame(table(seurat_subcluster_id))
}

P2_cell_number_df <- subcluster_cells_list[["P2"]]
P6_cell_number_df <- subcluster_cells_list[["P6"]]

# Rename the columns to apply merge
P2_cell_number_df <- P2_cell_number_df %>% dplyr::rename(Cluster = seurat_subcluster_id,
                                                  P2 = Freq)

P6_cell_number_df <- P6_cell_number_df %>% dplyr::rename(Cluster = seurat_subcluster_id,
                                                  P6 = Freq)

merged_p2_p6 <- merge(P2_cell_number_df, P6_cell_number_df, by = "Cluster", all = TRUE)

# Pivot to longer format for plotting
merged_p2_p6_long <- merged_p2_p6 %>%
  pivot_longer(cols = c(P6, P2), 
               names_to = "Sample", 
               values_to = "Cell_Count")

# P6 does not have any cells in subcluster 7 => replace NA with zero
merged_p2_p6_long <- merged_p2_p6_long %>% replace(is.na(.), 0)

# Add total number of cells per subcluster and their corresponding percentages
pct <- merged_p2_p6_long %>%
  group_by(Cluster) %>%
  mutate(Total_cells = sum(Cell_Count, na.rm = TRUE),
         Percentage = (Cell_Count / Total_cells) * 100) %>%
  ungroup()

# Visualize the distribution of the cells across the subclusters
ggplot(merged_p2_p6_long, aes(x = factor(Cluster), y = Cell_Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Cell_Count, y = Cell_Count + 100),
            color = "black", size = 4, position = position_dodge(width = 0.9)) +
  labs(x = "Clusters", y = "Cell Number") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 12),
    axis.title.y = element_text(size = 15, vjust = 0.9),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 0, hjust = 1)
  )


