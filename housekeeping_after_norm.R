### Check housekeeping profile after deseq control normalisation

counts_control_norm_long <- counts_control_norm %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_gene_id") %>%
  as_tibble() %>%
  pivot_longer(-ensembl_gene_id, names_to = "guide", values_to = "count") %>%
  inner_join(anno_hsap, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  inner_join(sample_metadata, by = c("guide" = "sample_name"))

housekeeping_counts_control_norm <- counts_control_norm_long %>%
  filter(external_gene_name %in% housekeeping_genes)

for (gene0 in housekeeping_genes){
  plot <- housekeeping_counts_control_norm %>%
    filter(external_gene_name == gene0) %>%
    ggplot(aes(guide, count)) +
    geom_point() +
    theme_cowplot(10) +
    theme(axis.text.x = element_text(angle = 45)) + 
    theme(axis.text.x = element_text(size = 4)) +
    ggtitle(gene0)
  assign(paste(gene0, "counts_plot", sep = '_'), envir = .GlobalEnv, plot)
}

plot_grid(ACTB_counts_plot, C1orf43_counts_plot, GAPDH_counts_plot, HMBS_counts_plot, PSMB4_counts_plot)

ggsave("housekeeping_counts_post_control_normalised.png", dpi = 300, width = 10, height = 5.4375, units = "in")

### DESeq2 auto normalised
counts_normal_norm_long <- counts_normal_norm %>%
  as.data.frame() %>%
  rownames_to_column("ensembl_gene_id") %>%
  as_tibble() %>%
  pivot_longer(-ensembl_gene_id, names_to = "guide", values_to = "count") %>%
  inner_join(anno_hsap, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  inner_join(sample_metadata, by = c("guide" = "sample_name"))

housekeeping_counts_normal_norm <- counts_normal_norm_long %>%
  filter(external_gene_name %in% housekeeping_genes)

for (gene0 in housekeeping_genes){
  plot <- housekeeping_counts_normal_norm %>%
    filter(external_gene_name == gene0) %>%
    ggplot(aes(guide, count)) +
    geom_point() +
    theme_cowplot(10) +
    theme(axis.text.x = element_text(angle = 45)) + 
    theme(axis.text.x = element_text(size = 4)) +
    ggtitle(gene0)
  assign(paste(gene0, "counts_plot", sep = '_'), envir = .GlobalEnv, plot)
}

plot_grid(ACTB_counts_plot, C1orf43_counts_plot, GAPDH_counts_plot, HMBS_counts_plot, PSMB4_counts_plot)

ggsave("housekeeping_counts_post_deseq2normal_normalised.png", dpi = 300, width = 10, height = 5.4375, units = "in")



### Ordered plot


pal_lancet()(5)
housekeeping_colours <- c("ACTB" = pal_lancet()(5)[1], 
                          "C1orf43" = pal_lancet()(5)[2], 
                          "GAPDH" = pal_lancet()(5)[3], 
                          "HMBS" = pal_lancet()(5)[4], 
                          "PSMB4" = pal_lancet()(5)[5])

housekeeping_counts_control_norm_sorted <- housekeeping_counts_control_norm %>%
  arrange(count)

housekeeping_counts_control_norm_centered_scaled_sorted <- housekeeping_counts_control_norm %>%
  group_by(external_gene_name) %>%
  mutate(count = count - mean(count)) %>%
  mutate(count = count / sd(count)) %>%
  arrange(count) 

all_housekeeping_control_norm <- housekeeping_counts_control_norm_centered_scaled_sorted %>%
  ggplot(aes(x = fct_reorder(guide, count), count, colour = external_gene_name)) + 
  geom_point() +
  scale_color_manual(values = housekeeping_colours) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 

housekeeping_noGapdh_control_norm <- housekeeping_counts_control_norm_centered_scaled_sorted %>%
  filter(external_gene_name != "GAPDH") %>%
  ggplot(aes(x = fct_reorder(guide, count), count, colour = external_gene_name)) + 
  geom_point() +
  scale_color_manual(values = housekeeping_colours) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 

housekeeping_noGapdh_noHmbs_control_norm <- housekeeping_counts_control_norm_centered_scaled_sorted %>%
  filter(external_gene_name != "GAPDH" & external_gene_name != "HMBS") %>%
  ggplot(aes(x = fct_reorder(guide, count), count, colour = external_gene_name)) + 
  geom_point() +
  scale_color_manual(values = housekeeping_colours) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

# plot_grid(all_housekeeping_control_norm, housekeeping_noGapdh_control_norm, housekeeping_noGapdh_noHmbs_control_norm)

all_housekeeping_control_norm

ggsave("all_housekeeping_postcontrol_normalisation.png", dpi = 300, width = 10, height = 5.4375, units = "in")

