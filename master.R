library(tidyverse)
library(ggrepel)
library(cowplot)
library(ggsci)

amplicon_counts <- read_delim(file = "amplicons_counts.csv", delim = ",") %>%
  mutate(gene = factor(gene))

amplicon_counts_long <- amplicon_counts %>%
  pivot_longer(cols = B2_S1:M9_S39, names_to = "guide", values_to = "count")

amplicon_counts_long %>%
  filter(gene == "C1orf43") %>%
  ggplot(aes(guide, count)) +
  geom_point() +
  theme_cowplot(10) +
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("C1orf43")

housekeeping_genes <- c("ACTB", "C1orf43", "HMBS", "PSMB4", "GAPDH")

housekeeping_counts <- amplicon_counts_long %>%
  filter(gene %in% housekeeping_genes)

for (gene0 in housekeeping_genes){
  plot <- housekeeping_counts %>%
    filter(gene == gene0) %>%
    ggplot(aes(guide, count)) +
    geom_point() +
    theme_cowplot(10) +
    theme(axis.text.x = element_text(angle = 45)) + 
    theme(axis.text.x = element_text(size = 4)) +
    ggtitle(gene0)
  assign(paste(gene0, "counts_plot", sep = '_'), envir = .GlobalEnv, plot)
}

plot_grid(ACTB_counts_plot, C1orf43_counts_plot, GAPDH_counts_plot, HMBS_counts_plot, PSMB4_counts_plot)

ggsave("housekeeping_counts.png", dpi = 300, width = 10, height = 5.4375, units = "in")

pal_lancet()(5)
housekeeping_colours <- c("ACTB" = pal_lancet()(5)[1], 
                          "C1orf43" = pal_lancet()(5)[2], 
                          "GAPDH" = pal_lancet()(5)[3], 
                          "HMBS" = pal_lancet()(5)[4], 
                          "PSMB4" = pal_lancet()(5)[5])

housekeeping_counts_sorted <- housekeeping_counts %>%
  arrange(count)

housekeeping_counts_centered_scaled_sorted <- housekeeping_counts %>%
  group_by(gene) %>%
  mutate(count = count - mean(count)) %>%
  mutate(count = count / sd(count)) %>%
  arrange(count) 

all_housekeeping <- housekeeping_counts_centered_scaled_sorted %>%
  ggplot(aes(x = fct_reorder(guide, count), count, colour = gene)) + 
  geom_point() +
  scale_color_manual(values = housekeeping_colours) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 

housekeeping_noGapdh <- housekeeping_counts_centered_scaled_sorted %>%
  filter(gene != "GAPDH") %>%
  ggplot(aes(x = fct_reorder(guide, count), count, colour = gene)) + 
  geom_point() +
  scale_color_manual(values = housekeeping_colours) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 

housekeeping_noGapdh_noHmbs <- housekeeping_counts_centered_scaled_sorted %>%
  filter(gene != "GAPDH" & gene != "HMBS") %>%
  ggplot(aes(x = fct_reorder(guide, count), count, colour = gene)) + 
  geom_point() +
  scale_color_manual(values = housekeeping_colours) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 

plot_grid(all_housekeeping, housekeeping_noGapdh, housekeeping_noGapdh_noHmbs)

ggsave("housekeeping_mean_var.png", dpi = 300, width = 10, height = 5.4375, units = "in")





summary_stats_per_gene <- amplicon_counts_long %>%
  group_by(gene) %>%
  summarise(mean = mean(count), var = var(count), n = n())

ggplot(summary_stats_per_gene, aes(mean, var, label = ifelse(mean > 4000, as.character(gene), ""))) + 
  geom_point() +
  theme_cowplot() + 
  geom_label_repel(label.size = 0.5) +
  coord_cartesian(xlim = c(0, 30000), ylim = c(0, 100000000))
  
amplicon_counts_long_log <- amplicon_counts_long %>%
  mutate(count = log2(count + 1))

summary_stats_per_gene_log <- amplicon_counts_long_log %>%
  group_by(gene) %>%
  summarise(mean = mean(count), var = var(count), n = n())

ggplot(summary_stats_per_gene_log, aes(mean, var, label = ifelse(var > 0.75 & mean > 4, as.character(gene), ""))) + 
  geom_point() +
  theme_cowplot() + 
  geom_label_repel(label.size = 0.5) + 
  ggtitle("Mean-Variance Relationship")


  