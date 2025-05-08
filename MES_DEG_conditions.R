##################################################################
### --------------- Step 6: DEG across conditions --------------- ###
##################################################################
print("running Step 6: DEG across conditions...")

#####################################################
### --------------- load packages --------------- ###
#####################################################
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(future)
library(tidyverse)
library(EnhancedVolcano)

#set working directory
setwd("..")

##################################################
### ---------- load data ---------- ###
##################################################

# load data after integration and PrepSCTFindMarkers()
dat <- readRDS("./rds/scanorama_res0.1.rds")

######################################################################
### ---------- identify DEG across conditions ---------- ###
######################################################################

plot_features <- function(dat, comparisons) {
  
  # Creating folders to house data
  main_dir <- getwd()
  table_dir <- "differential_expression"
  if(!dir.exists(file.path(main_dir, table_dir))) {
    dir.create(file.path(main_dir, table_dir))
  }

  # Loop through each comparison
  for (comparison in comparisons) {
    ident1 <- comparison$ident1
    ident2 <- comparison$ident2
    comp_label <- paste0(ident1, "_vs_", ident2)

    sub_dir <- file.path(table_dir, comp_label)
    if(!dir.exists(sub_dir)) {
      dir.create(sub_dir)
    }

    # By stimulation
    Idents(dat) <- dat$stim
    dat_stim <- subset(x=dat, idents = c(ident1, ident2))
  
    # Getting stim DEGs
    markers <- FindMarkers(dat_stim, ident.1 = ident1, 
                           ident.2 = ident2, 
                           min.pct = 0.2,
                           logfc.threshold = 0,
                           verbose = FALSE,
                           recorrect_umi = FALSE)
    write.csv(x = markers, file = paste0(sub_dir, "/stim_", comp_label, "_DEGs.csv"))

# plot volcano
# define the cutoffs for log fold change and adjusted P-value
logFCcutoff <- 0.1  # This value is typically set to 1 (log2 fold change), but you may adjust as needed.
pvalCutoff <- 0.05  # Adjusted P-value cutoff for significance.

# create custom key-value pairs for 'up', 'down', and 'ns' (not significant)
  keyvals <- ifelse(
    markers$avg_log2FC > logFCcutoff & markers$p_val_adj < pvalCutoff, 'red', #up
      ifelse(markers$avg_log2FC < -logFCcutoff & markers$p_val_adj < pvalCutoff, 'blue', #down
        'gray')) #ns
  names(keyvals)[keyvals == 'red'] <- 'up'
  names(keyvals)[keyvals == 'blue'] <- 'down'
  names(keyvals)[keyvals == 'gray'] <- 'ns'

# Filter significant markers based on cutoffs
significant_markers <- subset(markers, (abs(avg_log2FC) > logFCcutoff) & (p_val_adj < pvalCutoff))

# Function to safely select top markers
select_top_markers <- function(data, n = 5) {
  if (nrow(data) < n) {
    return(data) # Return all if less than n
  } else {
    return(head(data, n)) # Return top n
  }
}

# Sort and get top upregulated and downregulated markers
top_upregulated <- select_top_markers(significant_markers[order(-significant_markers$avg_log2FC), ])
top_downregulated <- select_top_markers(significant_markers[order(significant_markers$avg_log2FC), ])

# Combine top upregulated and downregulated
top_markers <- rbind(top_upregulated, top_downregulated)
  write.csv(x = top_markers, file = paste0(sub_dir, "/stim_", comp_label, "_top_DEGs.csv"))

# Create labels for top markers
labels <- rownames(top_markers)

# Create the volcano plot with custom colors and labels
volcano_plot <- EnhancedVolcano(markers,
                lab = rownames(markers),
                selectLab = labels,
                boxedLabels = TRUE,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = comp_label,
                pCutoff = pvalCutoff,
                FCcutoff = logFCcutoff,
                pointSize = 1.0,
                labSize = 3.0,
                colCustom = keyvals,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                max.overlaps = 50)
  ggsave(filename = paste0(sub_dir, "/stim_", comp_label, "_volcano.png"), plot = volcano_plot)
  
  # plot dot plot
  dotplot <- DotPlot(dat_stim, features = labels, cols = c("blue", "red")) + RotatedAxis()
  ggsave(filename = paste0(sub_dir, "/stim_", comp_label, "_dotplot.png"),
         width = 7, height = 5, plot = dotplot)

    # Loop through each cell type
    Idents(dat_stim) <- dat_stim$MES
    celltypes <- unique(dat_stim$MES)
    for (celltype in celltypes) {
      skip_to_next <- FALSE
      # Subset data for the current cell type
      dat_subset <- subset(x=dat_stim, idents = celltype)
      Idents(dat_subset) <- dat_subset$stim
      
      tryCatch({
        markers <- FindMarkers(dat_subset, ident.1 = ident1, 
                             ident.2 = ident2, 
                             min.pct = 0.2,
                             logfc.threshold = 0,
                             verbose = FALSE,
                             recorrect_umi = FALSE)
      }, error = function(e) {
        cat("Error: ", e$message, "\n")
        skip_to_next <<- TRUE
      })
      if(skip_to_next) { next }

      # Save DEGs
      write.csv(x = markers, file = paste0(sub_dir, "/", celltype, "_", comp_label, "_DEGs.csv"))

      # plot volcano
      # define the cutoffs for log fold change and adjusted P-value
      logFCcutoff <- 0.1  # This value is typically set to 1 (log2 fold change), but you may adjust as needed.
      pvalCutoff <- 0.05  # Adjusted P-value cutoff for significance.
      
      # create custom key-value pairs for 'up', 'down', and 'ns' (not significant)
      keyvals <- ifelse(
        markers$avg_log2FC > logFCcutoff & markers$p_val_adj < pvalCutoff, 'red', #up
        ifelse(markers$avg_log2FC < -logFCcutoff & markers$p_val_adj < pvalCutoff, 'blue', #down
               'gray')) #ns
      names(keyvals)[keyvals == 'red'] <- 'up'
      names(keyvals)[keyvals == 'blue'] <- 'down'
      names(keyvals)[keyvals == 'gray'] <- 'ns'
      
      # Filter significant markers based on cutoffs
      significant_markers <- subset(markers, (abs(avg_log2FC) > logFCcutoff) & (p_val_adj < pvalCutoff))
      
      # Function to safely select top markers
      select_top_markers <- function(data, n = 5) {
        if (nrow(data) < n) {
          return(data) # Return all if less than n
        } else {
          return(head(data, n)) # Return top n
        }
      }
      # Sort and get top upregulated and downregulated markers
      top_upregulated <- select_top_markers(significant_markers[order(-significant_markers$avg_log2FC), ])
      top_downregulated <- select_top_markers(significant_markers[order(significant_markers$avg_log2FC), ])
      # Combine top upregulated and downregulated
      top_markers <- rbind(top_upregulated, top_downregulated)
      # Create labels for top markers
      labels <- rownames(top_markers)

# Create the volcano plot with custom colors and labels
volcano_plot <- EnhancedVolcano(markers,
                lab = rownames(markers),
                selectLab = labels,
                boxedLabels = TRUE,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = paste0(celltype, "_", comp_label),
                pCutoff = pvalCutoff,
                FCcutoff = logFCcutoff,
                pointSize = 1.0,
                labSize = 3.0,
                colCustom = keyvals,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                max.overlaps = 50)
  ggsave(filename = paste0(sub_dir, "/", celltype, "_", comp_label, "_volcano.png"), plot = volcano_plot)
      
    }
  }

  return()
}

# Define comparisons
comparisons <- list(
  list(ident1 = "Sel", ident2 = "Control")
)

# Call the function
plot_features(dat, comparisons)

#####################################################
### --------------- Report run time --------------- ###
#####################################################

print("...finishing Step 6: DEG across conditions")