library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(knitr)
library(msigdbr)
library(fgsea)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(ComplexHeatmap)
library(colorRamp2)
library(tidyr)
library(patchwork)
library(ggtree)
library(gridExtra)
library(ape)
library(RColorBrewer)
library(reactome.db)
library(annotables)
library(biomaRt)
set.seed(5220)

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Define constants:

CONTEXT <- c("Sel vs Control")
INPUT_DIRS <- c(sprintf("/sb28_10x/sel-ctrl-tumor/differential_expression/gsea/csv"))
OUTPUT_DIR <- sprintf("/sb28_10x/sel-ctrl-tumor/differential_expression/gsea")

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Import the DESeq output:

data_list <- list()
for (dir in INPUT_DIRS) {
  file_names <- list.files(dir)
  for (i in file_names) {
    name <- gsub("_DEGs.csv", "", i)
    df <- read.table(paste(dir, i, sep = "/"), header = TRUE, sep = ",")
    colnames(df)[1] <- "feature"
    data_list[[name]] <- df
  }
}

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Prepare an ENSEMBL gene ID data list by converting using the `annotables` 
# package. Record the genes that end up going missing from ENSEMBL conversion.

data.list.symbols <- data_list

mapping.df <- annotables::grcm38
missing.genes <- list()
data.list.ensembl <- list()
for (name in names(data.list.symbols)) {
  deseq_table <- data.list.symbols[[name]]
  mapping.df <- mapping.df[!duplicated(mapping.df$symbol), ]
  mapping.table <- mapping.df[mapping.df$symbol %in% deseq_table$feature, ]
  mapping.table <- as.data.frame(mapping.table)
  rownames(mapping.table) <- mapping.table$symbol
  deseq_table$ensembl <- mapping.table[deseq_table$feature, "ensgene"]
  deseq_table.missing <- deseq_table[is.na(deseq_table$ensembl), ]
  missing.genes[[name]] <- deseq_table.missing$feature
  data.list.ensembl[[name]] <- deseq_table[!is.na(deseq_table$ensembl), ]
}
write.table(missing.genes[[1]], file = paste(OUTPUT_DIR, "missing_genes.txt", sep = "/"), row.names = FALSE, quote = FALSE)

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Remove any genes that don't pass DESeq's independent filtering, 
# indicated by having an NA p-value or adjusted p-value:

processed.deseq.tables <- list()
for (name in names(data.list.ensembl)) {
  deseq_table <- data.list.ensembl[[name]]
  deseq_table <- deseq_table[!is.na(deseq_table$p_val), ]
  deseq_table <- deseq_table[!is.na(deseq_table$p_val_adj), ]
  processed.deseq.tables[[name]] <- deseq_table
}
data.list <- processed.deseq.tables

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Generate an ordered list from the processed DESeq2 tables:

ranked_gene_set_list <- list()
for (name in names(processed.deseq.tables)) {
  deseq_table <- processed.deseq.tables[[name]]
  # Group by unique Ensembl ID and average
  deseq_table <- deseq_table %>%
    group_by(ensembl) %>%
    summarize(avg_log2FC = mean(avg_log2FC))
  deseq_table <- deseq_table %>% mutate(rank = rank(avg_log2FC,  ties.method = "random"))
  deseq_table <- deseq_table[order(-deseq_table$rank),] # Rank deterministically
  
  # Generate a ranked list
  gene_list <- deseq_table$avg_log2FC
  gene_names <- deseq_table$ensembl
  names(gene_list) <- gene_names
  ranked_gene_set_list[[name]] <- gene_list
}
## -----------------------------------------------------------------------------------------------------------------------------------------------
# Run GSEA using `msigdbr` and `fgsea`.

msigdbr_df <- msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)

set.seed(1) # required for deterministic performance
fgsea.output.list <- list()
for (name in names(ranked_gene_set_list)) {
  gsea_results <- fgsea(pathways = msigdbr_list, stats = ranked_gene_set_list[[name]],
                        maxSize = 500, eps = 0.0) # allow arbitrarily low p-values
  fgsea.output.list[[name]] <- gsea_results
}
saveRDS(fgsea.output.list, paste(OUTPUT_DIR, "fgsea_ensembl_ontology_list.rds", sep = "/"))

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Plot the enrichment results on a heatmap:

extracted_cols <- lapply(names(fgsea.output.list), function(name) {
  col <- fgsea.output.list[[name]]$NES
  names(col) <- fgsea.output.list[[name]]$pathway
  col
})
combined.df <- do.call(cbind, extracted_cols)
colnames(combined.df) <- names(fgsea.output.list)

ht = Heatmap(combined.df,
        name = "Enrichment Scores across stims",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_column_slices = FALSE,
        cluster_columns = TRUE,
        width = ncol(combined.df)*unit(4, "mm"),
        height = nrow(combined.df)*unit(4, "mm"),
)
pdf(paste(OUTPUT_DIR, "enrichment_combined_matrix.pdf", sep = "/"),
    width=10,
    height=15,
    useDingbats=FALSE)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Prepare the data for a bubble plot grouped within radiation conditions:

combined.df <- bind_rows(fgsea.output.list, .id = "MES")
combined.df$pathway <- gsub("HALLMARK_", "", combined.df$pathway)
combined.df$log10padj <- -log10(combined.df$padj)
unfiltered.df <- combined.df

combined.df <- combined.df %>%
  group_by(pathway) %>%
  filter(any(padj < 0.05)) %>%
  ungroup()

nes.df <- combined.df %>%
  dplyr::select(-pval, -padj, -log2err, -ES, -size, -leadingEdge, -log10padj)

nes.df <- nes.df %>%
  group_by(pathway, MES) %>%
  pivot_wider(names_from = MES, values_from = NES) %>%
  group_by(pathway)

nes.df <- as.data.frame(nes.df)
rownames(nes.df) <- nes.df$pathway
nes.df <- nes.df %>% dplyr::select(-pathway)

row_dist <- dist(nes.df, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")
row_order <- order.dendrogram(as.dendrogram(row_hclust))

df.ordered <- combined.df %>%
  arrange(match(pathway, rownames(nes.df)[row_order]), match(MES, colnames(nes.df)))

df.ordered$pathway <- factor(df.ordered$pathway, levels = rownames(nes.df)[row_order])
df.ordered$MES <- factor(df.ordered$MES, levels = colnames(nes.df))

# Add significance layer
df.ordered$Significance <- ifelse(df.ordered$padj < 0.05, "adjp < 0.05", "NA")

# Clean up names
df.ordered$Pathways <- df.ordered$pathway
levels(df.ordered$Pathways) <- levels(df.ordered$Pathways)

# Reorder Clusters based on the numeric part
#df.ordered$Clusters <- factor(df.ordered$Clusters, 
#                              levels = unique(df.ordered$Clusters[order(as.numeric(gsub("_Sel_vs_Control", "", df.ordered$Clusters)))]))


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Plot the bubble plot:

bubble.plot <- ggplot(df.ordered,
                      aes(x = MES)) +
  
  # Add the bubble and outline significance layers
  geom_point(aes(y = Pathways, size = log10padj, color = NES), alpha = 1.0) +
  geom_point(data = subset(df.ordered, padj < 0.05),
             aes(y = Pathways, size = log10padj, shape = 'adjp < 0.05'),
             color = "black", fill = NA, alpha = 0.8) +
  geom_point(data = subset(df.ordered, is.na(pval) | is.na(padj) | is.na(log10padj)), 
             aes(y = Pathways), 
             color = "black", size = 1, alpha = 0.6) +
  scale_size(range = c(2, 7)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey50"
  ) +
  scale_shape_manual(name = "Significance",
                     values = c(`adjp < 0.05` = 21),
                     guide = guide_legend(override.aes = list(color = "black", fill = NA))) +
  
  # Labels
  labs(
    color = "Effect size",
    size = "-Log10 padj"
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.height = unit(0.3, 'cm'),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank())

module.dendrogram <- ggtree(as.phylo(row_hclust)) + theme_tree()

vertical.plot <- bubble.plot
dendrogram.vertical.plot <- module.dendrogram
combined.plot <- dendrogram.vertical.plot | vertical.plot
combined.plot <- combined.plot + plot_layout(widths = c(0.2, 0.8)) + plot_annotation(
  title = "Enrichment of Hallmark pathways across all stims",
  subtitle = paste(CONTEXT, "Normalized")
)

ggsave(paste(OUTPUT_DIR, "enrichment_bubble_plot_fgsea_nes_within_all_stims.pdf", sep = "/"),
       combined.plot, height = 6, 
       width = 6, device = pdf)
## -----------------------------------------------------------------------------------------------------------------------------------------------
# Overlay mean LFCs on the existing plot

processed_deseq_tables <- list()
for (name in names(data.list.ensembl)) {
  deseq_table <- data.list.ensembl[[name]]
  deseq_table <- deseq_table[!is.na(deseq_table$p_val), ]
  deseq_table <- deseq_table[!is.na(deseq_table$p_val_adj), ]
  deseq_table <- deseq_table[deseq_table$p_val_adj < 0.05, ]
  deseq_table <- deseq_table[abs(deseq_table$avg_log2FC) > 0.1, ]
  processed_deseq_tables[[name]] <- deseq_table
}
data_list <- processed_deseq_tables

union_de_genes <- c()
for (name in names(data_list)) {
  deseq_table <- data_list[[name]]
  ensembl_ids <- deseq_table$ensembl
  union_de_genes <- unique(c(union_de_genes, ensembl_ids))
}

# Note that this list isn't exactly equal to the DE gene plots at the moment 
# because we also filter out all genes that don't have a valid ENSEMBL ID. 
# Next, we get the list of DESeq dataframes and filter them for these genes.
filter_to_de_genes <- function(deseq_table) {
  deseq_table <- subset(deseq_table, ensembl %in% union_de_genes)
}
de_gene_data_list <- lapply(data.list.ensembl, filter_to_de_genes)

# Generate an `msigdbr_df` that represents the genes representing each pathway.
msigdbr_df <- msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)

# For each perturbation x condition in the `de_gene_data_list`, 
# build a DESeq dataframe that includes only the intersection between the 
# pathway's genes and the current ones.
cluster_pathway_dfs <- list()
for (cluster in names(de_gene_data_list)) {
  deseq_df <- de_gene_data_list[[cluster]]
  for (pathway in names(msigdbr_list)) {
    pathway_intersected_df <- subset(deseq_df, ensembl %in% msigdbr_list[[pathway]])
    cluster_pathway <- paste(cluster, pathway, sep = "_")
    cluster_pathway_dfs[[cluster_pathway]] <- pathway_intersected_df
  }
}

# Find the average LFC across dataframe in `cluster_pathway_dfs`. 
# Note that for any gene that doesn't pass expression filters, we set its 
# contribution to the average to 0.
calculate_average <- function(deseq_table) {
  deseq_table$avg_log2FC[is.na(deseq_table$p_val) & is.na(deseq_table$p_val_adj)] <- 0
  return(mean(deseq_table$avg_log2FC))
}
cluster_pathway_means <- lapply(cluster_pathway_dfs, 
                                          calculate_average)

# Map the averages to the `df.ordered` table generated above.
df.ordered$cluster_pathway <- paste(df.ordered$MES, 
                                         df.ordered$Pathways, sep = "_")
df.ordered$cluster_pathway <- gsub(" ", "_", 
                                        df.ordered$cluster_pathway)

names(cluster_pathway_means) <- gsub(
  "HALLMARK_", 
  "", names(cluster_pathway_means))

df.ordered$cluster_pathway_means <- 
  unlist(cluster_pathway_means[df.ordered$cluster_pathway])

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Make the plot, but instead of using NES, use the mean:
bubble.plot <- ggplot(df.ordered,
                      aes(x = MES)) +
  
  # Add the bubble and outline significance layers
  geom_point(aes(y = Pathways, size = log10padj, color = cluster_pathway_means), alpha = 1.0) +
  geom_point(data = subset(df.ordered, padj < 0.05),
             aes(y = Pathways, size = log10padj, shape = 'adjp < 0.05'),
             color = "black", fill = NA, alpha = 0.8) +
  geom_point(data = subset(df.ordered, is.na(pval) | is.na(padj) | is.na(log10padj)), 
             aes(y = Pathways), 
             color = "black", size = 1, alpha = 0.6) +
  scale_size(range = c(2, 7)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey50",
    limits = c(-0.75, 0.75),
    oob = scales::squish
  ) +
  scale_shape_manual(name = "Significance",
                     values = c(`adjp < 0.05` = 21),
                     guide = guide_legend(override.aes = list(color = "black", fill = NA))) +
  
  # Labels
  labs(
    color = "Effect size",
    size = "-Log10 padj"
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.height = unit(0.3, 'cm'),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank())

module.dendrogram <- ggtree(as.phylo(row_hclust)) + theme_tree()

vertical.plot <- bubble.plot
dendrogram.vertical.plot <- module.dendrogram
combined.plot <- dendrogram.vertical.plot | vertical.plot
combined.plot <- combined.plot + plot_layout(widths = c(0.2, 0.8)) + plot_annotation(
  title = "Enrichment of Hallmark pathways across all stims",
  subtitle = paste(CONTEXT, "Normalized")
)
ggsave(paste(OUTPUT_DIR, "enrichment_bubble_plot_fgsea_nes_within_all_stims_logfc_overlay.pdf", sep = "/"),
       combined.plot, height = 6, 
       width = 6, device = pdf)
