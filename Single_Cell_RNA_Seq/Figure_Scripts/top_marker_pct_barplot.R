# Figure X: Analysis of top marker genes in late-passage cell populations and characterization of CRYAB interacting partners
# Barplot

# Input data: supplementary table 10

supplementary_table_10_df <- read.csv("/path/Supplementary Table 10.csv")

# Define cluster order based on their Passage status
subcluster_order <- c(0, 3, 4, 7, 1, 2, 5, 6)

# Define colors with the same palette as the other plots
ridge_colors <- c(
  "#bfadcc", "#ffe1f5", "#dfbbda", "#896790", 
  "#008080", "#006666", "#005959", "#004c4c"
)


supplementary_table_10_df %>%
  mutate(Clusters = factor(Clusters, levels = subcluster_order)) %>%
  ggplot(aes(x = Clusters, y = percentage, fill = Clusters)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = sprintf("%.1f", percentage)), vjust = -0.5) +
  # geom_text(aes(label = sprintf("%.1f%%", percentage)), vjust = -0.5) +
  labs(x = "Clusters", y = "Percentage (%)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  ) +
  scale_fill_manual(values = ridge_colors)
