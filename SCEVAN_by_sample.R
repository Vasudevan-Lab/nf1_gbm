####################################################
### --------------- load packages --------------- ###
#####################################################

library(spatstat.explore)
library(Seurat)
library(SCEVAN)
library(dplyr)

setwd("..")

# Creating folders to house data
main_dir <- getwd()
output_dir <- "SCEVAN_by_sample"
ifelse(!dir.exists(file.path(main_dir, output_dir)), 
       dir.create(file.path(main_dir, output_dir)), FALSE)

setwd("./SCEVAN_by_sample")

##################################################
### ---------- load data ---------- ###
##################################################
#load Seurat object
seurat_obj <- readRDS("../rds/NF1_patient_9_samples_integrated_annotated.rds")
DefaultAssay(seurat_obj)<-"RNA"
# add new column called tumor
seurat_obj$tumor <- ifelse(grepl("-like$", seurat_obj$celltype), "tumor", "normal")
# add new column called sample to rename orig.ident
orig_ident_to_sample <- c(
  "patient_1" = "S1",
  "patient_2" = "S9",
  "patient_3" = "S3",
  "patient_4" = "S2",
  "patient_5" = "S7",
  "patient_6" = "S4",
  "patient_7" = "S5",
  "patient_8" = "S6",
  "patient_10" = "S8"
)
seurat_obj@meta.data$sample <- orig_ident_to_sample[seurat_obj@meta.data$orig.ident]
Idents(seurat_obj) <- seurat_obj$sample
obj.list <- SplitObject(seurat_obj, split.by = "sample")

#####################################################
### --------------- functions --------------- ###
#####################################################
run_SCEVAN <- function(dat, sample_name) {
  # Creating folders to house data
  main_dir <- getwd()
  cnv_dir <- sample_name
  ifelse(!dir.exists(file.path(main_dir, cnv_dir)), 
         dir.create(file.path(main_dir, cnv_dir)), FALSE)
  setwd(cnv_dir)
  count_mtx <- GetAssayData(dat, slot = "counts")
  normal <- subset(dat, subset = tumor == "normal")
  norm_cells_vector <- colnames(normal)
  results <- SCEVAN::pipelineCNA(count_mtx, norm_cell = norm_cells_vector, sample = sample_name, par_cores = 50, SUBCLONES = TRUE, plotTree = TRUE)
  ## write to Seurat
  dat <- Seurat::AddMetaData(dat, metadata = results)
  setwd(main_dir)
  return(dat)
}

##################################################
### ---------- run SCEVAN and combine ---------- ###
##################################################
combined_obj <- lapply(names(obj.list), function(sample_name) {
  run_SCEVAN(obj.list[[sample_name]], sample_name = sample_name)
}) %>% 
  Reduce(function(x, y) merge(x, y), .)

saveRDS(combined_obj, file="combined_results.rds")
