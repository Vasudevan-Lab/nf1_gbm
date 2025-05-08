#########################################################################################################################################################################
#Cell type annotation
#scType https://github.com/IanevskiAleksandr/sc-type
library(scater)

#set working directory
main_dir <- getwd()

# load data
seurat_obj <- readRDS("../rds/scanorama_res0.15.rds")

# load libraries
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
#db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
db_ <- "../ref/ScTypeDB_full_mod_neftel.xlsx";
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
# can also include genes from cellMarker http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_download.html
gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = seurat_obj[["SCT"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
write.csv(sctype_scores, file = "../plots/scanorama/sctype_scores.csv")

seurat_obj@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_obj@meta.data$customclassif[seurat_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

seurat_obj[["sctype"]] = seurat_obj@meta.data$customclassif
umap_sctype<- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype')  
ggsave(filename = "../plots/scanorama/UMAP_sctype.png", plot = umap_sctype, width = 8, height = 5)
#########################################################################################################################################################################

# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL
write.csv(cL_resutls, file = "../plots/scanorama/cL_resutls.csv")
