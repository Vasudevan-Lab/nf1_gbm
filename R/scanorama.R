#####################################################
### --------------- Step 2: integration --------------- ###
#####################################################

# start Rscript timer
start_total_time <- Sys.time()
print("running Step 2: integraiton...")

#####################################################
### --------------- load packages --------------- ###
#####################################################
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(scCustomize)

library(reticulate)
python <- '/software/mambaforge/envs/scanpy/bin/python'

use_python(python)

#set working directory
setwd("..")

#DEBUG
# found that renamecells is outside of the if statement, so moved it inside

library(purrr)
Merge_Seurat_List <- function(
  list_seurat,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "SeuratProject"
) {
  # Check list_seurat is list
  if (!inherits(x = list_seurat, what = "list")) {
    cli_abort(message = "{.code list_seurat} must be environmental variable of class {.val list}")
  }

  # Check list_seurat is only composed of Seurat objects
  for (i in 1:length(x = list_seurat)) {
    if (!inherits(x = list_seurat[[i]], what = "Seurat")) {
      cli_abort("One or more of entries in {.code list_seurat} are not objects of class {.val Seurat}")
    }
  }

  # Check all barcodes are unique to begin with
  duplicated_barcodes <- list_seurat %>%
    lapply(colnames) %>%
    unlist() %>%
    duplicated() %>%
    any()

  if (duplicated_barcodes && is.null(x = add.cell.ids)) {
    cli_abort(message = c("There are overlapping cell barcodes present in the input objects",
                          "i" = "Please rename cells or provide prefixes to {.code add.cell.ids} parameter to make unique.")
    )
  }

  # Check right number of suffix/prefix ids are provided
  if (!is.null(x = add.cell.ids) && length(x = add.cell.ids) != length(x = list_seurat)) {
    cli_abort(message = "The number of prefixes in {.code add.cell.ids} must be equal to the number of objects supplied to {.code list_seurat}.")
    # Rename cells if provided
     list_seurat <- lapply(1:length(x = list_seurat), function(x) {
       list_seurat[[x]] <- RenameCells(object = list_seurat[[x]], add.cell.id = add.cell.ids[x])
     })
  }
  print(add.cell.ids[1])
 # list_seurat <- lapply(1:length(x = list_seurat), function(x) {
 #   list_seurat[[x]] <- RenameCells(object = list_seurat[[x]], add.cell.id = add.cell.ids[x])
 # })  

  # Merge objects
  merged_object <- reduce(list_seurat, function(x, y) {
    merge(x = x, y = y, merge.data = merge.data, project = project)
  })
}


##################################################
### ---------- load data ---------- ###
##################################################

#load data
dat <- readRDS("rds/merged_seurat_filtered.rds")

# Creating folders to house data
main_dir <- getwd()
plot_dir <- "plots"
ifelse(!dir.exists(file.path(main_dir, plot_dir)), 
       dir.create(file.path(main_dir, plot_dir)), FALSE)
sub_dir <- "scanorama"
ifelse(!dir.exists(file.path(main_dir, plot_dir, sub_dir)), 
       dir.create(file.path(main_dir, plot_dir, sub_dir)), FALSE)

######################################################################
### ---------- normalize and identify variable features ---------- ###
######################################################################

#Normalize each dataset individually
dat.list <-
  SplitObject(dat, split.by = "orig.ident")  %>%
  lapply(FUN = SCTransform, vst.flavor = "v2", return.only.var.genes = F, verbose = F)

assaylist <- list()
genelist <- list()
for(i in 1:length(dat.list))
{
  assaylist[[i]] <- t(as.matrix(GetAssayData(dat.list[[i]], "data")))
  genelist[[i]] <- rownames(dat.list[[i]])
}

scanorama <- import('scanorama')

# get sessioninfo, get Seurat and Scanorama versions
# sessionInfo() %>% print() %>% filter(grep("Seurat|scanorama", package, value = T)) %>% print()

integrated.data <- scanorama$integrate(assaylist, genelist)
corrected.data <- scanorama$correct(assaylist, genelist, return_dense=TRUE)
integrated.corrected.data <- scanorama$correct(assaylist, genelist, return_dimred=TRUE, return_dense=TRUE)

intdata <- lapply(integrated.corrected.data[[2]], t)
panorama <- do.call(cbind, intdata)
rownames(panorama) <- as.character(integrated.corrected.data[[3]])
colnames(panorama) <- unlist(sapply(assaylist, rownames))

intdimred <- do.call(rbind, integrated.corrected.data[[1]])
colnames(intdimred) <- paste0("PC_", 1:100)

#We also add standard deviations in order to draw Elbow Plots in Seurat

stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)

pan.seurat <- CreateSeuratObject(counts = panorama, assay = "pano",  project = "pano")

#Adding metadata from all previous objects 
pan.seurat@meta.data <- dat@meta.data

# VERY IMPORTANT: make sure that the rownames of your metadata slot 
# are the same as the colnames of your integrated expression matrix 

rownames(pan.seurat@meta.data) <- colnames(pan.seurat)
rownames(intdimred) <- colnames(pan.seurat)

#add RNA and SCT assays to combined.sct
combined.sct <- Merge_Seurat_List(dat.list)
# Get the count and data matrices for the pano assay
pano.counts <- GetAssayData(object = pan.seurat, assay = "pano", slot = "counts")
pano.data <- GetAssayData(object = pan.seurat, assay = "pano", slot = "data")
# Create a new assay object that combines the count and data matrices
pano.assay <- CreateAssayObject(counts = pano.counts)
pano.assay@data <- pano.data
# Add the new assay object to the Seurat object
combined.sct[["pano"]] <- pano.assay

DefaultAssay(combined.sct)='pano'
combined.sct[["pca"]] <- CreateDimReducObject(embeddings = intdimred, stdev = stdevs, key = "PC_", assay = "pano")
dims = 30
# Run a single integrated analysis on all cells
combined.sct <- FindNeighbors(combined.sct, dims = 1:dims, verbose = F)

# Cluster Tree
library(clustree)
dat<-combined.sct
resolution.range <- seq(from = 0, to = 0.5, by = 0.05)
dat <- FindClusters(object = dat, resolution = resolution.range, verbose = F)
cluster_tree<-clustree(dat, prefix = "pano_snn_res.")
ggsave(filename = "plots/scanorama/cluster_tree.png", width = 10, height = 10, plot = cluster_tree)

###########################################################
### ---------- save the integrated object ---------- ###
###########################################################

# set working directory to main folder
setwd(main_dir)

# save rds file
saveRDS(combined.sct, file="rds/scanorama.rds")

#####################################################
### --------------- Report run time --------------- ###
#####################################################

print("...finishing Step 2: integration")
#end timer
end_total_time <- Sys.time()
end_total_time - start_total_time
