library(tidyverse)
library(ggrepel)
library(cowplot)
library(ggsci)
library(DESeq2)
library(tximport)
library(ensembldb)
select <- dplyr::select
filter <- dplyr::filter
library(rhdf5)
library(biomaRt)

h5closeAll()

#Set up global parameters
skip_pcaExplorer = T
alpha = 0.01
res_lfc = 0 #lfc threshold for DESeq2 results functions
lfc = 0 #lfc threshold for post-results filtering genes
replicate_replace = Inf
dpi = 300
width = 7.25 #Dimensions set to span entire width of A4 publication, 4:3 ratio
height = 5.4375
units = "in"

# kallisto_path <- "/zfs/analysis/SWH_CRISPRi/191120-2_targeted_rna_seq_batch_1_/fastq/kallisto_20200131/"
kallisto_path_v1 <- "/zfs/analysis/SWH_CRISPRi/202003_updated_GAPDH/kallisto_amplicon_v1"

#updated GAPDH reference
kallisto_path_v2 <- "/zfs/analysis/SWH_CRISPRi/202003_updated_GAPDH/kallisto_amplicon_v2"
tx2gene <- read_delim("202003_updated_GAPDH/tx2gene", delim = " ")

# txdb <- makeTxDbFromGFF('Homo_sapiens.GRCh38.99.gtf', organism = "Homo sapiens")
# tx2gene <- AnnotationDbi::select(txdb, keys(txdb, keytype = "TXNAME"), "GENEID", "TXNAME")

extra_info <- read_delim("sample_index_swh_crispri.csv", delim =",")

sample_metadata <- read_delim("/zfs/analysis/SWH_CRISPRi/191120-2_targeted_rna_seq_batch_1_/list.txt", delim = "\t", col_names = "sample_name") %>%
  mutate(sample_code = sample_name) %>%
  mutate(sample_name = str_extract(sample_name, "[^_]+")) %>%
  inner_join(extra_info, by = c("sample_name" = "sample_name")) %>%
  mutate(gene_target = recode(gene_target, PARK2 = "PRKN")) %>%
  mutate(batch = ifelse(gene_target == "LRRK2", "batch_1", 
                        ifelse(gene_target == "PINK1", "batch_1", 
                        ifelse(gene_target == "PARK2", "batch_1",
                        ifelse(gene_target == "PARK7", "batch_1",
                        ifelse(gene_target == "SNCA", "batch_1",
                        ifelse(gene_target == "CHCHD2", "batch_1", "batch_2")))))))

#Housekeeping definitions
all_housekeeping_genes <- c("ACTB", "C1orf43", "HMBS", "PSMB4", "GAPDH")

good_housekeeping <- c("ACTB", "C1orf43", "PSMB4")
good_housekeeping <- c("ACTB", "C1orf43", "PSMB4", "GAPDH")


###Examine original counts file
# amplicon_counts <- read_delim(file = "amplicons_counts.csv", delim = ",") %>%
#   mutate(gene = factor(gene))
# 
# amplicon_counts_long <- amplicon_counts %>%
#   pivot_longer(cols = B2_S1:M9_S39, names_to = "guide", values_to = "count")
# 
# amplicon_counts_long %>%
#   filter(gene == "C1orf43") %>%
#   ggplot(aes(guide, count)) +
#   geom_point() +
#   theme_cowplot(10) +
#   theme(axis.text.x = element_text(angle = 90)) + 
#   ggtitle("C1orf43")
# 
# amplicon_counts_long %>%
#   filter(gene == "PINK1") %>%
#   ggplot(aes(guide, count)) +
#   geom_point() +
#   theme_cowplot(10) +
#   theme(axis.text.x = element_text(angle = 90)) + 
#   ggtitle("PINK1")

###Summary stats from original counts file

# summary_stats_per_gene <- amplicon_counts_long %>%
#   group_by(gene) %>%
#   summarise(mean = mean(count), var = var(count), n = n())
# 
# ggplot(summary_stats_per_gene, aes(mean, var, label = ifelse(mean > 4000, as.character(gene), ""))) + 
#   geom_point() +
#   theme_cowplot() + 
#   geom_label_repel(label.size = 0.5) +
#   coord_cartesian(xlim = c(0, 30000), ylim = c(0, 100000000))
#   
# amplicon_counts_long_log <- amplicon_counts_long %>%
#   mutate(count = log2(count + 1))
# 
# summary_stats_per_gene_log <- amplicon_counts_long_log %>%
#   group_by(gene) %>%
#   summarise(mean = mean(count), var = var(count), n = n())
# 
# ggplot(summary_stats_per_gene_log, aes(mean, var, label = ifelse(var > 0.75 & mean > 4, as.character(gene), ""))) + 
#   geom_point() +
#   theme_cowplot() + 
#   geom_label_repel(label.size = 0.5) + 
#   ggtitle("Mean-Variance Relationship")


### Start using kallisto counts

##prepare data for pca on all samples
name = 'swh_batch_1_all'
design = ~ 1
min_counts = 3

#create list of kallisto input files
files <- file.path(kallisto_path_v2, sample_metadata$sample_code, "abundance.h5")
names(files) <- sample_metadata$sample_name
#import kallisto files
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = T)
#create deseq2 object
dds <- DESeqDataSetFromTximport(txi, colData = sample_metadata, design = design)
#filter low counts
keep_feature <- rowMeans(counts(dds)) >= min_counts
dds_removed <- dds[!keep_feature, ]
dds <- dds[keep_feature, ]

# keep_housekeeping <- rownames(dds) %in% good_housekeeping$ensembl_gene_id
# dds_housekeeping <- dds[keep_housekeeping, ]
# dds_housekeeping <- DESeq(dds_housekeeping, minReplicatesForReplace = replicate_replace)
# 
# estimateSizeFactors(dds_housekeeping)
# mean(colMeans(normalizationFactors(dds_housekeeping)))

which(rownames(dds) %in% good_housekeeping)
dds_control_norm <- estimateSizeFactors(dds, controlGenes = which(rownames(dds) %in% good_housekeeping))
mean(colMeans(normalizationFactors(dds_control_norm)))
range(colMeans(normalizationFactors(dds_control_norm)))
dds_control_norm <-DESeq(dds_control_norm, minReplicatesForReplace = replicate_replace)
counts_control_norm <- counts(dds_control_norm, normalized = T) 
vsd_control_norm <- varianceStabilizingTransformation(dds_control_norm, blind = T)
vsd_control_norm_long <- assay(vsd_control_norm) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  pivot_longer(-gene, names_to = "guide", values_to = "count") %>% 
  inner_join(sample_metadata, by = c("guide" = "sample_name"))


counts_no_norm <- counts(dds, normalized = F)

vsd_no_norm <- varianceStabilizingTransformation(dds, blind = T)
vsd_no_norm_long <- assay(vsd_no_norm) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  pivot_longer(-gene, names_to = "guide", values_to = "count") %>% 
  inner_join(sample_metadata, by = c("guide" = "sample_name"))

dds_normal_norm <- estimateSizeFactors(dds)
mean(colMeans(normalizationFactors(dds_normal_norm)))
range(colMeans(normalizationFactors(dds_normal_norm)))
dds_normal_norm <-DESeq(dds_normal_norm, minReplicatesForReplace = replicate_replace)
counts_normal_norm <- counts(dds_normal_norm, normalized = T) 

###PCA

plot_pca <- function(pca_df, group, label, title) {
  if (label == T) {
    group <- enquo(group)
    var_explained <- pca_df$sdev^2/sum(pca_df$sdev^2) #Calculate PC variance
    pca_df$x %>% 
      as.data.frame() %>%
      rownames_to_column("sample_name") %>%
      inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
      ggplot(aes(x=PC1,y=PC2, label = !! group, color = !! group)) + 
      geom_point(aes(color = !! group), size = 4) +
      labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
           y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
      geom_label_repel(box.padding = 0.5) +
      theme_bw() +
      # scale_color_lancet() + 
      # scale_fill_lancet() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
      ggtitle(title) +
      theme(legend.position = "none")
  } else {
    group <- enquo(group)
    var_explained <- pca_df$sdev^2/sum(pca_df$sdev^2) #Calculate PC variance
    pca_df$x %>% 
      as.data.frame() %>%
      rownames_to_column("sample_name") %>%
      inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
      ggplot(aes(x=PC1,y=PC2, color = !! group)) + 
      geom_point(aes(color = !! group), size = 4) +
      labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
           y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
      theme_bw() +
      # scale_color_lancet() + 
      # scale_fill_lancet() + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) + 
      ggtitle(title) +
      theme(legend.position = "none")
  }
}

#PCA all IP samples

#No normalisation
pca_no_norm <- prcomp(t(assay(vsd_no_norm)), scale=T) #Calculate PCs
var_explained_no_norm <- pca_no_norm$sdev^2/sum(pca_no_norm$sdev^2) #Calculate PC variance
plot_pca(pca_no_norm, group = batch, label = T, title = 'All samples')

#label no normalisation by GAPDH count
GAPDH_counts <- vsd_no_norm_long %>% 
  filter(gene == "GAPDH") %>%
  select(guide, count)
sample_metadata <- sample_metadata %>%
  inner_join(GAPDH_counts, by = c("sample_name" = "guide")) %>% 
  rename("count" = "gapdh_vsd_count")
pca_no_norm$x %>% 
  as.data.frame() %>%
  rownames_to_column("sample_name") %>%
  inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
  ggplot(aes(x=PC1,y=PC2, label = gene_target, color = gapdh_vsd_count)) + 
  geom_point(aes(color = gapdh_count), size = 4) +
  labs(x=paste0("PC1: ",round(var_explained_no_norm[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained_no_norm[2]*100,1),"%")) +
  geom_label_repel(box.padding = 0.5) +
  theme_cowplot() +
  # scale_color_lancet() + 
  # scale_fill_lancet() +
  ggtitle("Non-normalised PC1 and PC2 biplot") +
  labs(color = "GAPDH Count")
  
ggsave("housekeeping_counts_prenormalisation_kallisto.png", dpi = 300, width = 10, height = 5.4375, units = "in")

#Control normalised
pca_control_norm <- prcomp(t(assay(vsd_control_norm)), scale=T) #Calculate PCs
var_explained_control_norm <- pca_control_norm$sdev^2/sum(pca_control_norm$sdev^2) #Calculate PC variance
plot_pca(pca_control_norm, group = batch, label = T, title = 'All samples')

ggsave(filename = "pca_control_norm.png", width = width*2, height = height*2, dpi = dpi, units = units) 

#DESeq2 auto normalised
pca_normal_norm <- prcomp(counts_normal_norm, scale=T) #Calculate PCs
var_explained_normal_norm <- pca_normal_norm$sdev^2/sum(pca_normal_norm$sdev^2) #Calculate PC variance

plot_pca(pca_normal_norm, group = batch, label = T, title = 'All samples')

ggsave(filename = "pca_normal_norm.png", width = width, height = height, dpi = dpi, units = units) 

###Try removing PC1
counts_control_norm_batch_removed <- limma::removeBatchEffect(counts_control_norm, sample_metadata$batch)

rv <- rowVars(counts_control_norm_batch_removed) #Calculate row variance
rv_order <- order(rv, decreasing = TRUE)
keep <- head(rv_order, max(1, nrow(counts_control_norm_batch_removed)*(1)))
matrix_high_var <- counts_control_norm_batch_removed[keep, ] %>% #Select and transpose top n vsd genes
  t()
pca_control_norm_batch_removed <- prcomp(matrix_high_var, scale=T) #Calculate PCs
var_explained_control_norm_batch_removed <- pca_control_norm_batch_removed$sdev^2/sum(pca_control_norm_batch_removed$sdev^2) #Calculate PC variance

plot_pca(pca_control_norm_batch_removed, group = gene_target, label = T, title = 'All samples')
ggsave(filename = "pca_control_norm_batch_corrected.png", width = width*2, height = height*2, dpi = dpi, units = "in") 


###Get external gene names

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")

anno_hsap <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'ensembl_gene_id', values = rownames(counts_control_norm), mart = mart)

# 
# ###Create control normalised counts table
# 
# counts_control_norm %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   rename(ensembl_gene_id = "") %>%
#   # inner_join(anno_hsap, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
#   pivot_longer(-ensembl_gene_id, names_to = "guide", values_to = "count") 

###Housekeeping counts based on original counts

# housekeeping_counts_original <- amplicon_counts_long %>%
#   filter(gene %in% all_housekeeping_genes)
# 
# for (gene0 in all_housekeeping_genes){
#   plot <- housekeeping_counts_original %>%
#     filter(gene == gene0) %>%
#     ggplot(aes(guide, count)) +
#     geom_point() +
#     theme_cowplot(10) +
#     theme(axis.text.x = element_text(angle = 45)) + 
#     theme(axis.text.x = element_text(size = 4)) +
#     ggtitle(gene0)
#   assign(paste(gene0, "counts_plot", sep = '_'), envir = .GlobalEnv, plot)
# }
# 
# plot_grid(ACTB_counts_plot, C1orf43_counts_plot, GAPDH_counts_plot, HMBS_counts_plot, PSMB4_counts_plot)
# 
# ggsave("housekeeping_counts_prenormalisation.png", dpi = 300, width = 10, height = 5.4375, units = "in")

###Housekeeping counts based on kallisto counts

counts_no_norm_long <- counts_no_norm %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  pivot_longer(-gene, names_to = "guide", values_to = "count") %>% 
  inner_join(sample_metadata, by = c("guide" = "sample_name"))

counts_control_norm_long <- counts_control_norm %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  pivot_longer(-gene, names_to = "guide", values_to = "count") %>%
  inner_join(sample_metadata, by = c("guide" = "sample_name"))

housekeeping_counts_no_norm <- counts_no_norm_long %>%
  filter(gene %in% all_housekeeping_genes)

for (gene0 in all_housekeeping_genes){
  plot <- housekeeping_counts_no_norm %>%
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

ggsave("housekeeping_counts_prenormalisation_kallisto.png", dpi = 300, width = 10, height = 5.4375, units = "in")

### Control normalised counts
housekeeping_counts_control_norm <- counts_control_norm_long %>%
  filter(gene %in% all_housekeeping_genes)

for (gene0 in all_housekeeping_genes){
  plot <- housekeeping_counts_control_norm %>%
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

# ggsave("housekeeping_counts_prenormalisation_kallisto.png", dpi = 300, width = 10, height = 5.4375, units = "in")

###Sorted housekeeping original counts
pal_lancet()(5)
housekeeping_colours <- c("ACTB" = pal_lancet()(5)[1], 
                          "C1orf43" = pal_lancet()(5)[2], 
                          "GAPDH" = pal_lancet()(5)[3], 
                          "HMBS" = pal_lancet()(5)[4], 
                          "PSMB4" = pal_lancet()(5)[5])

housekeeping_counts_no_norm_sorted <- housekeeping_counts_no_norm %>%
  arrange(count)

housekeeping_counts_no_norm_centered_scaled_sorted <- housekeeping_counts_no_norm %>%
  group_by(gene) %>%
  mutate(count = count - mean(count)) %>%
  mutate(count = count / sd(count)) %>%
  arrange(count) 

all_housekeeping <- housekeeping_counts_no_norm_centered_scaled_sorted %>%
  ggplot(aes(x = fct_reorder(guide, count), count, colour = gene)) + 
  geom_point() +
  scale_color_manual(values = housekeeping_colours) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 

housekeeping_noGapdh <- housekeeping_counts_no_norm_centered_scaled_sorted %>%
  filter(gene != "GAPDH") %>%
  ggplot(aes(x = fct_reorder(guide, count), count, colour = gene)) + 
  geom_point() +
  scale_color_manual(values = housekeeping_colours) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 

housekeeping_noGapdh_noHmbs <- housekeeping_counts_no_norm_centered_scaled_sorted %>%
  filter(gene != "GAPDH" & gene != "HMBS") %>%
  ggplot(aes(x = fct_reorder(guide, count), count, colour = gene)) + 
  geom_point() +
  scale_color_manual(values = housekeeping_colours) +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

plot_grid(all_housekeeping, housekeeping_noGapdh, housekeeping_noGapdh_noHmbs)

ggsave("housekeeping_prenormalisation_sorted.png", dpi = 300, width = 10, height = 5.4375, units = "in")

all_housekeeping

ggsave("all_housekeeping_prenormalisation.png", dpi = 300, width = 10, height = 5.4375, units = "in")

###Sorted housekeeping kallisto counts

housekeeping_counts_control_norm_sorted <- housekeeping_counts_control_norm %>%
  arrange(count)

housekeeping_counts_control_norm_centered_scaled_sorted <- housekeeping_counts_control_norm_sorted %>%
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

plot_grid(all_housekeeping_control_norm, housekeeping_noGapdh_control_norm, housekeeping_noGapdh_noHmbs_control_norm)

ggsave("housekeeping_prenormalisation_sorted_kallisto.png", dpi = 300, width = 10, height = 5.4375, units = "in")

all_housekeeping_control_norm

ggsave("all_housekeeping_prenormalisation_kallisto.png", dpi = 300, width = 10, height = 5.4375, units = "in")

####Gene-level plots
genes_targeted <- counts_control_norm_long %>%
  select(gene_target) %>%
  arrange() %>%
  distinct()

genes_measured <- counts_control_norm_long %>%
  select(external_gene_name) %>%
  arrange() %>%
  distinct() %>%
  as.list()

counts_control_norm_long %>% 
  select(gene_target) %>%
  arrange() %>%
  distinct()

medians <- counts_control_norm_long %>%
  filter(gene != gene_target) %>%
  group_by(gene) %>%
  summarise(median = median(count)) %>%
  mutate(log2median = log2(median + 1))



gene_targeted = "PARK7"

counts_control_norm_long %>%
  inner_join(medians, by = c("gene" = "gene")) %>%
  filter(gene_target == gene_targeted) %>%
  mutate(log2count = (log2(count) + 1)) %>%
  mutate(is_gene_targeted = ifelse(gene == gene_target, 'y', 'n')) %>%
  select(guide, 
         count, 
         log2count, 
         gene, 
         gene_target, 
         median, 
         log2median, 
         is_gene_targeted) %>%
  pivot_longer(cols = c("log2count", "log2median"), names_to = "count_type", values_to = "log_value") %>% 
  ggplot(aes(x = gene, y = log_value, colour = count_type)) + 
  geom_point() + 
  theme_cowplot() + 
  scale_colour_npg() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 6)) +
  # theme(legend.position = "None") + 
  ggtitle(gene_targeted)


for(target in genes_targeted$gene_target){
  plot <- counts_control_norm_long %>%
    filter(gene == target) %>%
    mutate(log2count = (log2(count) + 1)) %>%
    mutate(is_gene_targeted = ifelse(gene_target == target, 'y', 'n')) %>%
    ggplot(aes(x = gene, log2count, colour = is_gene_targeted)) +
    # geom_violin(aes(fill = is_gene_targeted)) + 
    geom_point(position = position_jitterdodge()) +
    theme_cowplot() +
    theme(legend.position = "None") + 
    scale_color_lancet() +
    scale_fill_lancet() +
    theme(axis.title.x = element_blank())
  assign(paste(target, "counts_plot", sep = '_'), envir = .GlobalEnv, plot)
}

plot_grid(LRRK2_counts_plot, GBA_counts_plot, GRN_counts_plot, PINK1_counts_plot, PRKN_counts_plot, VPS35_counts_plot, PARK7_counts_plot, ATP13A2_counts_plot, SNCA_counts_plot, ATP6V0A1_counts_plot)

ggsave(filename = "target_genes_scatter.png", width = width*2, height = height*2, dpi = dpi, units = "in")


#Per gene measured plots
for(gene in genes_measured$external_gene_name){
  plot <- counts_control_norm_long %>%
    filter(external_gene_name == gene) %>%
    mutate(log2count = (log2(count + 1))) %>%
    mutate(is_gene_targeted = ifelse(gene_target == gene, 'y', 'n')) %>%
    ggplot(aes(x = external_gene_name, log2count, colour = gene_target)) +
    # geom_violin(aes(fill = is_gene_targeted)) + 
    geom_jitter() +
    theme_cowplot() +
    theme(legend.position = "None") 
    # scale_color_lancet() +
    # scale_fill_lancet()
  assign(paste(gene, "measured_counts_plot", sep = '_'), envir = .GlobalEnv, plot)
}

PRKN_measured_counts_plot


for(gene in genes_targeted$gene_target){
  diff <- counts_control_norm_long %>% 
  filter(gene_target == gene & external_gene_name == gene) %>%
  mutate(log2count = log2(count + 1)) %>%
  inner_join(medians, by = c("external_gene_name" = "external_gene_name")) %>%
  mutate(diff_to_median = log2count - log2median) %>%
    ggplot(aes(x = external_gene_name, diff_to_median, colour = gene)) +
    geom_jitter() + 
    theme_cowplot() + 
    theme(legend.position = "None")
  assign(paste(gene, "diff_to_median_plot", sep = '_'), envir = .GlobalEnv, diff)
  
           }

plot_grid(LRRK2_diff_to_median_plot, 
          GBA_diff_to_median_plot, 
          GRN_diff_to_median_plot, 
          PINK1_diff_to_median_plot, 
          PRKN_diff_to_median_plot, 
          VPS35_diff_to_median_plot, 
          PARK7_diff_to_median_plot, 
          ATP13A2_diff_to_median_plot, 
          SNCA_diff_to_median_plot, 
          ATP6V0A1_diff_to_median_plot)

individual_diff <- counts_control_norm_long %>% 
  select(ensembl_gene_id, count, external_gene_name, gene_target) %>%
  filter(gene_target == external_gene_name) %>%
  mutate(log2count = log2(count + 1)) %>%
  inner_join(medians, by = c("external_gene_name" = "external_gene_name")) %>%
  mutate(diff_to_median = log2count - log2median) %>%
  mutate(diff_proportion = 2^diff_to_median)

  counts_control_norm_long %>% 
  select(ensembl_gene_id, count, external_gene_name, gene_target) %>%
  filter(gene_target == external_gene_name) %>%
  mutate(log2count = log2(count + 1)) %>%
  inner_join(medians, by = c("external_gene_name" = "external_gene_name")) %>%
  mutate(diff_to_median = log2count - log2median) %>%
    group_by(gene_target) %>%
    summarise(mean_diff = mean(diff_to_median), sd_diff = sd(diff_to_median)) %>%
    mutate(mean_diff_proportion = 2^mean_diff, 
           lower_CI_proportion = 2^(mean_diff - (sd_diff * 1.96)), 
           upper_CI_proportion  = 2^(mean_diff + (sd_diff * 1.96))) %>%
    arrange(mean_diff_proportion) %>%
    ggplot(aes(x = fct_reorder(gene_target, mean_diff_proportion), y = mean_diff_proportion, fill = gene_target)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = lower_CI_proportion, ymax = upper_CI_proportion), width=.2,
                    position=position_dodge(.9)) + 
    theme_cowplot(14) + 
    scale_fill_npg() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
    theme(axis.title.x = element_blank()) + 
    theme(legend.position = "None") +
    ylab("Proportion of median expression") +
    geom_text(
      aes(label = paste(round((mean_diff_proportion*100), digits = 2), "%", sep = "")),
      position = position_dodge(0.9),
      vjust = -15,
      face = "bold",
      size = 3
    ) +
    geom_point(data = individual_diff, aes(x = gene_target, y = diff_proportion))
  
  ggsave(filename = "target_genes_barplot_with_points.png", width = width, height = height, dpi = dpi, units = "in")

  

heatmap(counts_control_norm, scale = "row")

library("gplots")
heatmap.2(counts_control_norm, scale = "row", col = bluered(100), 
          trace = "none", density.info = "none")

library("pheatmap")
pheatmap(counts_control_norm, cutree_cols = 8)

library(ComplexHeatmap)

rownames(counts_control_norm)

match(colnames(counts_control_norm), sample_metadata$sample_name)

test <- counts_control_norm %>% 
  t() %>% 
  scale() %>% 
  t()

rownames(test) <- anno_hsap$external_gene_name
rownames(test)
colnames(test) <- sample_metadata$gene_target
colnames(test)

Heatmap(test, show_row_names = T)
