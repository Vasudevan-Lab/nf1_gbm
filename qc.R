#####################################################
### --------------- Step 1: QC --------------- ###
#####################################################

#start Rscript timer
start_total_time <- Sys.time()
print("running Step 1: QC...")


#####################################################
### --------------- load packages --------------- ###
#####################################################
library(Seurat)
library(magrittr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(future)

#set working directory
setwd("..")

##################################################
### ---------- load data ---------- ###
##################################################

read_data <- function(path_to_data, file_name, sample_name) {
  dat <- Read10X(paste0(path_to_data, "/", file_name, "/outs/filtered_feature_bc_matrix"))
  object <- CreateSeuratObject(counts = dat, 
                               min.cells = 3,
                               min.features = 200,
                               project = sample_name)
  return(object)
}

path_to_folder <- '/raleighlab/data1/mirchiak'

# Excluded patients
excluded <- c(9, 11, 12, 13, 14)

# Create a list to store Seurat objects
patients_list <- list()

# Loop through patient numbers
for (i in 1:14) {
  if (!(i %in% excluded)) {
    sample_name <- paste("patient", i, sep = "_")
    file_name <- paste("KM_NF1_10x", i, sep = "_")
    patients_list[[sample_name]] <- read_data(path_to_data = path_to_folder, file_name = file_name, sample_name = sample_name)
  }
}

#directory to house qc data
main_dir <- getwd()
qc_dir <- "qc"
ifelse(!dir.exists(file.path(main_dir, qc_dir)), 
       dir.create(file.path(main_dir, qc_dir)), FALSE)

# Merge Seurat objects
# First, make sure the list isn't empty (this can cause errors when merging)
if (length(patients_list) == 0) {
    stop("No Seurat objects found in the list!")
}
merged_seurat <- merge(x = patients_list[[1]], y = patients_list[2:length(patients_list)], add.cell.ids = names(patients_list))

# species assignment (species <- "human" or species <- "mouse")
species <- "human" 

# patterns based on species
if (species == "human") {
  mt_pattern <- "MT-"
  rb_pattern <- "RP[SL]"
} else if (species == "mouse") {
  mt_pattern <- "mt-"
  rb_pattern <- "rp[sl]"
} else {
  stop("Unknown species type!")
}

merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = mt_pattern)
merged_seurat[["percent.rb"]] <- PercentageFeatureSet(merged_seurat, pattern = rb_pattern)

# pre_filtering violin plots
filter = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
qc_plots <- VlnPlot(merged_seurat, features = filter, ncol = 2)
ggsave(filename = "qc/qc_plots_pre_filter.png", plot = qc_plots, width = 20, height = 10)

###############################################
### ---------- data qc filtering ---------- ###
###############################################
qc <- function(dat, name) {
  
  #start timer
  start_time <- Sys.time()
  print("running qc filtering...")
  
  #parallel
  plan("multisession", workers = 20)
  
  # Creating folders to house data
  main_dir <- getwd()
  sub_dir <- name
  qc_dir <- "qc"
  ifelse(!dir.exists(file.path(main_dir, qc_dir)), 
         dir.create(file.path(main_dir, qc_dir)), FALSE)
  ifelse(!dir.exists(file.path(main_dir, qc_dir, sub_dir)), 
         dir.create(file.path(main_dir, qc_dir, sub_dir)), FALSE)
  
  # Filter the data
  
  dat[["percent.mt"]] <-
    PercentageFeatureSet(dat, pattern = mt_pattern)
  dat[["percent.rb"]] <-
    PercentageFeatureSet(dat, pattern = rb_pattern)
  
  filter = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
  qc_plots <- VlnPlot(dat, features = filter, ncol = 2)
  ggsave(filename = paste0("qc/", name, "/qc_plots.png"), plot = qc_plots)
  
  filter_table <- list()
  filter_value = c("min", "mean", "median", "max",  "sd" , "mad")
  dat.metadata <- dat@meta.data[, c(filter)]
  
  for (i in 1:length(filter_value)) {
    filter_table[[i]] = sapply(dat.metadata, filter_value[i])
  }
  
  filter_table <- data.frame(filter_table)
  colnames(filter_table) <- filter_value
  write.csv(x = filter_table, file = paste0("qc/", name, "/filter_table_",
                                            name, ".csv"))
  
  dat <-
    subset(
      x = dat,
      subset = nFeature_RNA > 500 &
        nFeature_RNA < as.integer(filter_table["nFeature_RNA", "median"] + 3 * filter_table["nFeature_RNA", "mad"]) &
        nCount_RNA > 1000 &
        nCount_RNA < as.integer(filter_table["nCount_RNA", "median"] + 3 * filter_table["nCount_RNA", "mad"]) &
        percent.mt < 10 &
        percent.rb < (filter_table["percent.rb", "median"] + 3 * filter_table["percent.rb", "mad"]),
      return.null = TRUE)
  
  if (is.null(dat)) {
    warning("No cells found for after filtering!")
  } else {
    qc_filtered_plots <- VlnPlot(dat, features = filter, ncol = 2)
    ggsave(filename = paste0("qc/", name, "/qc_filtered_plots.png"), plot = qc_filtered_plots)
  }
  
  #end timer
  end_time <- Sys.time()
  end_time - start_time
  
  return(dat)
}

# Create a list to store filtered Seurat objects
patients_list_filtered <- list()

# Loop through patient numbers
for (i in 1:14) {
  if (!(i %in% excluded)) {
    sample_name <- paste("patient", i, sep = "_")
    patients_list_filtered[[sample_name]] <- qc(patients_list[[sample_name]], sample_name)
  }
}

# Merge Seurat objects after filtering
# First, make sure the list isn't empty (this can cause errors when merging)
if (length(patients_list_filtered) == 0) {
    stop("No Seurat objects found in the list!")
}
merged_seurat_filtered <- merge(x = patients_list_filtered[[1]], y = patients_list_filtered[2:length(patients_list_filtered)], add.cell.ids = names(patients_list_filtered))

merged_seurat_filtered[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_filtered, pattern = mt_pattern)
merged_seurat_filtered[["percent.rb"]] <- PercentageFeatureSet(merged_seurat_filtered, pattern = rb_pattern)

# post_filtering violin plots
filter = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")
qc_plots <- VlnPlot(merged_seurat_filtered, features = filter, ncol = 2)
ggsave(filename = "qc/qc_plots_post_filter.png", plot = qc_plots, width = 20, height = 10)

#cell count table
count_table <- rbind(table(merged_seurat_filtered@meta.data$orig.ident))
rownames(count_table) <- c("n")
write.csv(x = count_table, file = "qc/post_filtering_count_table.csv")

#####################################################
### --------------- save data --------------- ###
#####################################################

# Creating folders to house data
main_dir <- getwd()
data_dir <- "rds"
ifelse(!dir.exists(file.path(main_dir, data_dir)), 
       dir.create(file.path(main_dir, data_dir)), FALSE)

saveRDS(merged_seurat_filtered, file="rds/merged_seurat_filtered.rds")

#####################################################
### --------------- Report run time --------------- ###
#####################################################

print("...finishing Step 1: QC")
#end timer
end_total_time <- Sys.time()
end_total_time - start_total_time