# Boxplot of quantitative fluorescence intensity data
# Date: 16-05-2025

beta_gal_data <- read_excel("path/to/betagal_updated_2025.xlsx")

beta_gal_data <- beta_gal_data %>% dplyr::rename("Passage 2" = "P2",
                                       "Passage 6" = "P6")

beta_gal_longer <- beta_gal_data %>% pivot_longer(cols = c(`Passage 2`, `Passage 6`), names_to = "scores")

ggplot(beta_gal_longer, aes(x = scores, y = value, fill = scores)) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.7) +
  stat_boxplot(geom = "errorbar",
               width = 0.15) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5, color = "black") +
  theme_classic() +
  labs(title = "",
       y = "Normalized Fluorescence Intensity \n (Arbitrary Units)",
       x = "") +
  scale_fill_manual(values = c("#45087b", "#00ff00")) +
  theme(
    axis.text = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.title.x = element_text(size = 25),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(size = 30),
    plot.title = element_text(size = 30),
    legend.position = "none")
