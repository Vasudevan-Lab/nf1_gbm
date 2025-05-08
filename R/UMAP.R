#####################################################
### --------------- Step 4: UMAP --------------- ###
#####################################################

# start Rscript timer
start_total_time <- Sys.time()
print("running Step 4: UMAP...")

#####################################################
### --------------- load packages --------------- ###
#####################################################
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(reticulate)
library(scCustomize)

#set working directory
setwd("..")
main_dir <- getwd()

##################################################
### ---------- load data ---------- ###
##################################################

#load data
dat <- readRDS("./rds/scanorama.rds")

########################################
### ---------- Clustering ---------- ###
########################################
resolution = 0.15
# Run a single integrated analysis on all cells
dat <- FindClusters(dat, resolution = resolution, verbose = F) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE)

######################################################
### ---------- Visualization with UMAP ---------- ###
######################################################

# set working directory to plot directory
plot_dir<- "./plots/scanorama"
setwd(paste0(main_dir, "/", plot_dir))

library(RColorBrewer)
num_clusters <- length(unique(dat$seurat_clusters))
if (num_clusters <= 12) {
    ClusterColors <- brewer.pal(12, "Paired")
} else {
    ClusterColors <- colorRampPalette(brewer.pal(12, "Paired"))(num_clusters)
}

# UMAP by cluster
umap_clusters <- DimPlot(dat, reduction = "umap", group.by = "seurat_clusters", cols =  ClusterColors, label = T)
ggsave(filename = "UMAP_seurat_clusters.png", plot = umap_clusters, width = 7, height = 7)

# UMAP by sample
umap_samples <- DimPlot(dat, reduction = "umap", group.by="orig.ident")
ggsave(filename = "UMAP_samples.png", plot = umap_samples, width = 7, height = 7)

# Get number of cells per cluster and per sample of origin
count_cell_per_cluster <- table(dat@meta.data$orig.ident, dat@meta.data$seurat_clusters)
write.csv(x = count_cell_per_cluster, file = "count_per_sample_per_cluster.csv")

####################################################
### ---------- plot selected features ---------- ###
####################################################

# Features
DefaultAssay(dat)='SCT'

# UMAP cell cycle
library(stringr)
s <- str_to_title(cc.genes.updated.2019$s.genes)
g2.m <- str_to_title(cc.genes.updated.2019$g2m.genes)
dat <- CellCycleScoring(dat,
                        g2m.features = g2.m,
                        s.features = s, 
                        set.ident = TRUE)
umap_cell_cycle <- DimPlot(dat, reduction = "umap", group.by = "Phase")
ggsave(filename = "UMAP_cell_cycle.png", plot = umap_cell_cycle, width = 7, height = 7)

# Feature plot to check bias
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
qc_features_plot <- FeaturePlot(dat, features = qc_features, ncol = 2)
ggsave(filename = "qc_features_plot.png", plot = qc_features_plot, width = 7, height = 7)

###########################################################
### ---------- CD45 ---------- ###
###########################################################
# Creating folders to house data
feature_dir <- "./plots/scanorama/features"
ifelse(!dir.exists(file.path(main_dir, feature_dir)), 
       dir.create(file.path(main_dir, feature_dir)), FALSE)
setwd(paste0(main_dir, "/", feature_dir))

# plot PTPRC
features <- "PTPRC"
DefaultAssay(dat)='SCT'
Idents(dat) <- dat$seurat_clusters 

feature_plot <- FeaturePlot(dat, features = features, label = T)
ggsave(filename = "CD45_feature_plot.png", plot = feature_plot, width = 7, height = 7)

violin_plot_clusters <- VlnPlot(dat, features = features, cols = ClusterColors)
ggsave(filename = "CD45_violin_plot_clusters.png", plot = violin_plot_clusters, width = 20, height = 5)

###########################################################
### ---------- silhouette score ---------- ###
###########################################################
# set working directory to plot directory
setwd(paste0(main_dir, "/", plot_dir))

# Calculate silhouette score

library(cluster)
reduction <- "pca"
dims <- 1:30
Idents(dat)<-dat$seurat_clusters
dist.matrix <- dist(x = Embeddings(object = dat[[reduction]])[, dims])
clusters <- dat$seurat_clusters
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)

# prepare colors:
clust.col <- ClusterColors
sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
# plot:

png(filename = "silhouette.png", width = 5, height = 25, units="in", res = 200)
par(cex=0.8)
plot(sil,
     main = paste(length(unique(clusters)), "clusters"),
     border = sil.cols,
     col = sil.cols,
     do.col.sort = FALSE) 
while (!is.null(dev.list()))
  dev.off()

##################################################
### ---------- scMRMA annotation ---------- ###
##################################################
library(scMRMA)

## Run scMRMA annotation
DefaultAssay(dat)='RNA'
Idents(dat)<-dat$seurat_clusters
result <- scMRMA(input=dat, species="Hs",selfClusters=Idents(dat))
# species= "Hs" (default) or "Mm"
# selfClusters Use fixed clusters in each level

## UMAP plot with scMRMA annotation
dat[["scMRMA"]] <- result$multiR$annotationResult[colnames(dat),ncol(result$multiR$annotationResult)]
umap_scMRMA <- DimPlot(dat,reduction = "umap",group.by = "scMRMA",label = TRUE,repel = TRUE)
ggsave(filename = "UMAP_scMRMA.png", width = 7, height = 5, plot = umap_scMRMA)

# save a dataframe of scMRMA result for each cluster

file <- cbind(result$uniformR$annotationResult, result$uniformR$meta)
#delete rownames of file
rownames(file) <- NULL
# shorten file by deleting repeating rows and sort by cluster
file <- file %>%
  distinct() %>%
  mutate(clusters = as.numeric(gsub("Root_", "", clusters)) - 1) %>%
  arrange(clusters)

write.csv(file, file = "scMRMA_result.csv")

######################################################################
### ---------- identify cell type markers per cluster ---------- ###
######################################################################

# set working directory to main folder
setwd(main_dir)

# Creating folders to house data
table_dir <- "differential_expression"
ifelse(!dir.exists(file.path(main_dir, table_dir)), 
       dir.create(file.path(main_dir, table_dir)), FALSE)

print("running PrepSCTFindMarkers and FindAllMarkers for cluster markers...")

# run PrepSCTFindMarkers to ensure fixed value is set properly
Idents(dat)<-'seurat_clusters'
DefaultAssay(dat)='SCT'
# debug: Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) :
options(future.globals.maxSize = 8000 * 1024^2)
dat <- PrepSCTFindMarkers(dat, verbose = F)

# set working directory to differential expression folder
setwd(paste0(main_dir, "/differential_expression"))

# Getting DEGs for each cluster vs everything else
markers_clust <- FindAllMarkers(dat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, min.diff.pct = 0.2)
markers_clust %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(x = markers_clust, file = "clusters_DEGs.csv")

top10 <- markers_clust %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(x = top10, file = "clusters_top10_DEGs.csv")

# Heatmap of cluster markers
heatmap_clusters <- DoHeatmap(dat, features = top10$gene, group.colors = ClusterColors) + NoLegend()
ggsave(filename = "heatmap_clusters.png", plot = heatmap_clusters, width = 25, height = 15)

###########################################################
### ---------- save the integrated object ---------- ###
###########################################################

# set working directory to main folder
setwd(main_dir)

# save integrated object as rds file
saveRDS(dat, file="rds/scanorama_res0.15.rds")

#####################################################
### --------------- Report run time --------------- ###
#####################################################

print("...finishing Step 4: UMAP")
#end timer
end_total_time <- Sys.time()
end_total_time - start_total_time