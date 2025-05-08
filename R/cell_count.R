#####################################################
### --------------- cell count analysis --------------- ###
#####################################################

# start Rscript timer
start_total_time <- Sys.time()
print("running Step: cell count analysis...")

#####################################################
### --------------- load packages --------------- ###
#####################################################
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)

#set working directory
setwd("..")

##################################################
### ---------- load data ---------- ###
##################################################

# load data after integration and PrepSCTFindMarkers()
dat <- readRDS("./rds/scanorama_res0.15.rds")

# Creating folders to house data
main_dir <- getwd()
stats_dir <- "./plots/scanorama/counts_stats"
ifelse(!dir.exists(file.path(main_dir, stats_dir)), 
       dir.create(file.path(main_dir, stats_dir)), FALSE)

##############################################
### ---------- plot cell percentage ---------- ###
##############################################

setwd(stats_dir)

library(RColorBrewer)
num_clusters <- length(unique(dat$seurat_clusters))
if (num_clusters <= 12) {
    ClusterColors <- brewer.pal(12, "Paired")
} else {
    ClusterColors <- colorRampPalette(brewer.pal(12, "Paired"))(num_clusters)
}
#counts plot by sample
sample_counts_plot<- ggplot(dat@meta.data, aes(x=orig.ident, fill=seurat_clusters)) + 
  geom_bar(color="black", position = "fill", width = 0.8) +
  labs(fill='Cell Type') +
  scale_fill_manual(values = ClusterColors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        legend.text = element_text(color="black",size=13),legend.title = element_text(color="black",size=13),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15), axis.title=element_text(size=15), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+ 
  labs(y ="Composition of cells", x= NULL)+
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)
ggsave(filename = "counts_sample.png", plot = sample_counts_plot, width = 7, height = 5)

#counts plot by cluster
cluster_counts_plot<- ggplot(dat@meta.data, aes(x=seurat_clusters, fill=orig.ident)) + 
  geom_bar(color="black", position = "fill", width = 0.8) +
  labs(fill='Cell Type') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        legend.text = element_text(color="black",size=13),legend.title = element_text(color="black",size=13),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15), axis.title=element_text(size=15), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+ 
  labs(y ="Composition of cells", x= NULL)+
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)
ggsave(filename = "counts_cluster.png", plot = cluster_counts_plot, width = 7, height = 5)

celltype_counts_plots <- ggplot(dat@meta.data, aes(seurat_clusters, fill=seurat_clusters))+
  geom_bar(stat="count",colour = "black",width = 0.8)+  
  scale_fill_manual(values = ClusterColors)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        axis.line = element_line(colour = "black"), axis.text.y = element_blank(),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=10),
        axis.text=element_text(size=6.5), axis.title=element_text(size=15),
        plot.margin = unit(c(-0.1, 0.5,2.5, -0.5),"cm"))+labs(y ="number of cells", x= NULL)+ 
  scale_y_continuous(expand=c(0,0),trans ="log2")+coord_flip()+
  scale_x_discrete(limits=rev)+
  theme(legend.position = "none")
ggsave(filename = "counts_cell_per_celltype.png", plot = celltype_counts_plots, width = 7, height = 7)

cluster_counts_plot<- ggplot(dat@meta.data, aes(x=seurat_clusters, fill=Phase)) + 
  geom_bar(color="black", position = "fill", width = 0.8) +
  labs(fill='Cell Type') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        legend.text = element_text(color="black",size=13),legend.title = element_text(color="black",size=13),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15), axis.title=element_text(size=15), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+ 
  labs(y ="Composition of cells", x= NULL)+
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)
ggsave(filename = "counts_phase_cluster.png", plot = cluster_counts_plot, width = 7, height = 5)

# SCEVAN
cluster_counts_plot<- ggplot(dat@meta.data, aes(x=seurat_clusters, fill=class)) + 
  geom_bar(color="black", position = "fill", width = 0.8) +
  labs(fill='Cell Type') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        legend.text = element_text(color="black",size=13),legend.title = element_text(color="black",size=13),
        axis.line = element_line(colour = "black"), axis.text.y = element_text(color="black",size=12),axis.text.x = element_text(color="black",angle = 45,hjust = 1,size=15),
        axis.text=element_text(size=15), axis.title=element_text(size=15), 
        plot.margin = unit(c(0.2, 0.5,0.2, 0.5),"cm"))+ 
  labs(y ="Composition of cells", x= NULL)+
  scale_y_continuous(expand = c(0,0),limits = c(0,1),breaks = c(0,0.20,0.40,0.60,0.80,1),labels = scales::percent)
ggsave(filename = "counts_tumor_cluster.png", plot = cluster_counts_plot, width = 7, height = 5)