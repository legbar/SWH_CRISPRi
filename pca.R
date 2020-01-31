
##prepare data for pca on all samples
name = 'swh_batch_1_all'
design = ~ 1
min_counts = 3

#create list of kallisto input files
files <- file.path(kallisto_path, sample_metadata$sample_code, "abundance.h5")
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

which(rownames(dds) %in% good_housekeeping$ensembl_gene_id)
dds_control_norm <- estimateSizeFactors(dds, controlGenes = c(7, 36, 43))
mean(colMeans(normalizationFactors(dds_control_norm)))
range(colMeans(normalizationFactors(dds_control_norm)))
dds_control_norm <-DESeq(dds_control_norm, minReplicatesForReplace = replicate_replace)
counts_control_norm <- counts(dds_control_norm, normalized = T) 


dds_normal_norm <- estimateSizeFactors(dds)
mean(colMeans(normalizationFactors(dds_normal_norm)))
range(colMeans(normalizationFactors(dds_normal_norm)))
dds_normal_norm <-DESeq(dds_normal_norm, minReplicatesForReplace = replicate_replace)
counts_normal_norm <- counts(dds_normal_norm, normalized = T) 



#PCA

plot_pca <- function(pca_df, group, label, title) {
  if (label == T) {
    group <- enquo(group)
    var_explained <- pca_df$sdev^2/sum(pca_df$sdev^2) #Calculate PC variance
    pca_df$x %>% 
      as.data.frame() %>%
      rownames_to_column("sample_name") %>%
      inner_join(sample_metadata, by = c("sample_name" = "sample_name")) %>%
      ggplot(aes(x=PC1,y=PC2, label = sample_name, color = !! group)) + 
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

#Control normalised
rv <- rowVars(counts_control_norm) #Calculate row variance
rv_order <- order(rv, decreasing = TRUE)
keep <- head(rv_order, max(1, nrow(counts_control_norm)*(1)))
matrix_high_var <- counts_control_norm[keep, ] %>% #Select and transpose top n vsd genes
  t()
pca_control_norm <- prcomp(matrix_high_var, scale=T) #Calculate PCs
var_explained_control_norm <- pca$sdev^2/sum(pca$sdev^2) #Calculate PC variance

plot_pca(pca_control_norm, group = guide_group, label = T, title = 'All samples')

ggsave(filename = "pca_control_norm.png", width = width, height = height, dpi = dpi, units = units) 

#DESeq2 auto normalised
rv <- rowVars(counts_normal_norm) #Calculate row variance
rv_order <- order(rv, decreasing = TRUE)
keep <- head(rv_order, max(1, nrow(counts_normal_norm)*(1)))
matrix_high_var <- counts_normal_norm[keep, ] %>% #Select and transpose top n vsd genes
  t()
pca_normal_norm <- prcomp(matrix_high_var, scale=T) #Calculate PCs
var_explained_normal_norm <- pca$sdev^2/sum(pca$sdev^2) #Calculate PC variance

plot_pca(pca_normal_norm, group = pca_group, label = T, title = 'All samples')

ggsave(filename = "pca_normal_norm.png", width = width, height = height, dpi = dpi, units = units) 

###

counts_control_norm %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  rename(ensembl_gene_id = rowname) %>%
  # inner_join(anno_hsap, by = c("ensembl_gene_id" = "ensembl_gene_id")) %>%
  pivot_longer(-ensembl_gene_id, names_to = "guide", values_to = "count") 

