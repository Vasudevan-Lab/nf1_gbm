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

#set working directory
setwd("..")
main_dir <- getwd()

##################################################
### ---------- load data ---------- ###
##################################################

#load data
dat <- readRDS("./rds/NF1_patient_9_samples_integrated_annotated.rds")

##################################################
### ---------- curated signatures ---------- ###
##################################################
# for 2, 8, 9, 10, 11, 13

olig_markers <- c("PLP1","MBP")
mg_markers <- c("CD14", "CD68", "ITGAM", "CD163", "ITGA4", "ITGAX")
neuron_markers <- c("SYT1","RBFOX3","GRIK1","ATP1B1")
tcell_markers <- c("PTPRC", "CD2","IL7R", "ITGA4")
endothelial_markers <- c("CD34","VWF")

markers <- c(olig_markers, mg_markers, neuron_markers, tcell_markers, endothelial_markers)

setwd("./plots/features")
violin_plot <- VlnPlot(dat,
        features = markers,
        split.by = "seurat_clusters",
        cols = ClusterColors,
        flip = TRUE,
        stack = TRUE)
ggsave(filename="curated_markers_violin_plot.png", plot = violin_plot, width = 7, height = 7)


##################################################
### ---------- myeloid signatures from 2017 paper ---------- ###
##################################################
# load xls file from 2017 paper
signatures <- readxl::read_excel("./code-paper-v2/13059_2017_1362_MOESM5_ESM.xlsx")
MG_markers <- na.omit(signatures$MG_Markers)
Mac_markers <- na.omit(signatures$Mac_Markers)

setwd("./plots/features") 

feature_plot <- FeaturePlot(dat, features = MG_markers, label = T)
ggsave(filename="2017_MG_feature_plot.png", plot = feature_plot, width = 10, height = 10)

feature_plot <- FeaturePlot(dat, features = Mac_markers, label = T)
ggsave(filename="2017_Mac_feature_plot.png", plot = feature_plot, width = 10, height = 40)

library(RColorBrewer)
num_clusters <- length(unique(dat$seurat_clusters))
if (num_clusters <= 12) {
    ClusterColors <- brewer.pal(12, "Paired")
} else {
    ClusterColors <- colorRampPalette(brewer.pal(12, "Paired"))(num_clusters)
}

# Pre-define the indexes for clusters with specific color themes
MES <- c(1, 5)
AC <- 3
OPC <- 0
NPC <- 6

# Directly assign specified colors to the identified clusters
ClusterColors[c(MES + 1, AC + 1, OPC + 1, NPC + 1)] <- c("#A6CEE3", "#3385BB", # Light and darker blue for MES
                                                         "#6DBD57",             # Green for AC
                                                         "#FE8D19",             # Orange for OPC
                                                         "#E42622")             # Red for NPC

# Identify the rest of the clusters that need their colors randomized
rest_of_clusters <- setdiff(0:(num_clusters - 1), c(MES, AC, OPC, NPC))

# Randomize the rest of the colors
set.seed(123) 
ClusterColors[rest_of_clusters + 1] <- sample(ClusterColors[rest_of_clusters + 1])

# Generate the violin plot
violin_plot <- VlnPlot(dat, 
        features = MG_markers,
        split.by = "seurat_clusters",
        cols = ClusterColors,
        flip = TRUE,
        stack = TRUE)
ggsave(filename="2017_MG_violin_plot.png", plot = violin_plot, width = 10, height = 10)

violin_plot <- VlnPlot(dat, 
        features = Mac_markers,
        split.by = "seurat_clusters",
        cols = ClusterColors,
        flip = TRUE,
        stack = TRUE)
ggsave(filename="2017_Mac_violin_plot.png", plot = violin_plot, width = 10, height = 40)

microglia_genes <- c("NAV3", "P2RY12", "CX3CR1", "TMEM119", "SERPINE1", "SYNDIG1", "TAL1", "SALL1", "WIPF3")
bmdm_genes <- c("FPR3", "VEGFA", "TGFBI", "ITGA4", "VCAN", "TNS1")

violin_plot <- VlnPlot(dat,
        features = microglia_genes,
        split.by = "seurat_clusters",
        cols = ClusterColors,
        flip = TRUE,
        stack = TRUE)
ggsave(filename="microglia_violin_plot.png", plot = violin_plot, width = 10, height = 10)

violin_plot <- VlnPlot(dat,
        features = bmdm_genes,
        split.by = "seurat_clusters",
        cols = ClusterColors,
        flip = TRUE,
        stack = TRUE)
ggsave(filename="bmdm_violin_plot.png", plot = violin_plot, width = 10, height = 5)

#module score feature plot
gene_lists <- list(
  MG = MG_markers,
  Mac = Mac_markers,
  BMDM = bmdm_genes)

# Loop over the gene lists
for (name in names(gene_lists)) {
  # Add Module Score
  dat <- AddModuleScore(dat, features = list(gene_lists[[name]]), name = name)

  # violin plot for clusters
    violin_plot_clusters <- VlnPlot(dat, features = paste0(name, "1"), cols = ClusterColors, pt.size = 0)
    ggsave(filename = paste0(name, "_violin_plot_clusters.png"), plot = violin_plot_clusters, width = 20, height = 5)
}

# stacked violin
violin_plot <- VlnPlot(dat, features = c("MG1", "Mac1", "BMDM1"), 
                       split.by = "seurat_clusters", cols = ClusterColors,
        flip = TRUE,
        stack = TRUE)
ggsave(filename="stacked_violin_plot.png", plot = violin_plot, width = 10, height = 5)

##################################################
### ---------- myeloid signatures ---------- ###
##################################################

setwd("./plots/features") 

microglia_genes <- c("NAV3", "P2RY12", "CX3CR1", "TMEM119", "SERPINE1", "SYNDIG1", "TAL1", "SALL1", "WIPF3")
bmdm_genes <- c("FPR3", "VEGFA", "TGFBI", "ITGA4", "VCAN", "TNS1")

feature_plot <- FeaturePlot(dat, features = microglia_genes, label = T)
ggsave(filename="microglia_feature_plot.png", plot = feature_plot, width = 10, height = 10)

feature_plot <- FeaturePlot(dat, features = bmdm_genes, label = T)
ggsave(filename="bmdm_feature_plot.png", plot = feature_plot, width = 10, height = 10)

# Generate the violin plot
violin_plot <- VlnPlot(dat, features = microglia_genes, pt.size = 0) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename="microglia_violin_plot.png", plot = violin_plot, width = 10, height = 10)

##################################################
### ---------- dendritic signatures ---------- ###
##################################################
dc_genes <- c("ITGAX", "CD86", "CD80", "ITGAM", "AIF1", "HLA-DQB1")
violin_plot <- VlnPlot(dat,
        features = dc_genes,
        split.by = "seurat_clusters",
        cols = ClusterColors,
        flip = TRUE,
        stack = TRUE)
ggsave(filename="dc_violin_plot.png", plot = violin_plot, width = 10, height = 5)


##########################################
### ---------- data by primary/recurrent ---------- ###
##########################################

#create meta data by condition based on Oncoplot_v2.xlsx
dat2 <- dat@meta.data$orig.ident
dat2 <- replace(dat2, dat2=="patient_1", "Primary")
dat2 <- replace(dat2, dat2=="patient_2", "Recurrent")
dat2 <- replace(dat2, dat2=="patient_3", "Recurrent")
dat2 <- replace(dat2, dat2=="patient_4", "Primary")
dat2 <- replace(dat2, dat2=="patient_5", "Recurrent")
dat2 <- replace(dat2, dat2=="patient_6", "Recurrent")
dat2 <- replace(dat2, dat2=="patient_7", "Primary")
dat2 <- replace(dat2, dat2=="patient_8", "Primary")
dat2 <- replace(dat2, dat2=="patient_10", "Primary")
dat@meta.data$sample_type <- dat2

##########################################
### ---------- data by celltype ---------- ###
##########################################

# Create a named vector for the cluster-to-cell type mapping
cluster_to_celltype <- c(
  "0" = "OPC-like", 
  "1" = "MES-like", 
  "2" = "Oligodendrocytes", 
  "3" = "AC-like", 
  "4" = "Microglia", 
  "5" = "MES-like", 
  "6" = "NPC-like", 
  "7" = "Microglia", 
  "8" = "Neurons", 
  "9" = "Neurons", 
  "10" = "Neurons", 
  "11" = "T cells", 
  "12" = "Microglia", 
  "13" = "Endothelia"
)

# Extract the cluster IDs
clusters <- Idents(dat)

# Rename the clusters using the mapping
celltypes <- cluster_to_celltype[as.character(clusters)]
dat$celltype <- celltypes

###########################################################
### ---------- save the integrated object ---------- ###
###########################################################

# set working directory to main folder
setwd(main_dir)

# save integrated object as rds file
saveRDS(dat, file="rds/NF1_patient_9_samples_integrated_annotated.rds")