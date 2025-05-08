library(infercnv)
library(Seurat)

setwd("..")
# Creating folders to house data
main_dir <- getwd()
sub_dir <- "cnv"
ifelse(!dir.exists(file.path(main_dir, sub_dir)), 
       dir.create(file.path(main_dir, sub_dir)), FALSE)

seurat_obj <- readRDS("rds/scanorama_res0.15.rds")
# subset to exclude CD45+ clusters 4,7,11,12
seurat_obj <- subset(seurat_obj, idents = c(0,1,2,3,5,6,8,9,10,13))
DefaultAssay(seurat_obj)<-"RNA"
Idents(seurat_obj) <- seurat_obj$seurat_clusters
setwd("./cnv")

counts_matrix = GetAssayData(seurat_obj, slot="counts")
ann <- cbind(colnames(seurat_obj), as.character(seurat_obj@meta.data$seurat_clusters))
write.table(ann, file = "annotation_file_endo.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix,
                                    annotations_file="annotation_file_endo.txt",
                                    delim="\t",
                                    gene_order_file="ref/infercnv/human_gene_ordering_file.txt",
                                    ref_group_names=c("13")) 

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="infercnv_endo", 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             num_threads = 50,
                             HMM=TRUE)

seurat_obj = infercnv::add_to_seurat(seurat_obj,
                                     infercnv_output_path="infercnv_endo",
                                     top_n=10)

saveRDS(seurat_obj, file="infercnv_endo.rds")