library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readxl)
library(nichenetr)
library(msigdbr)

setwd("..")
main_dir<-getwd()
#load data 
dat <- readRDS("./rds/tumor.rds")

dat$MES <- ifelse(dat$Quadrant == "MESlike", "MES", "nonMES")

##########################################
### ---------- genesets ---------- ###
##########################################
feature_dir <- "neftel/MES/genesets_more"
ifelse(!dir.exists(file.path(main_dir, feature_dir)), 
       dir.create(file.path(main_dir, feature_dir)), FALSE)
setwd(file.path(main_dir, feature_dir))

# Helper function to perform Wilcoxon test and return results
perform_wilcox_test <- function(dat, formula) {
    # Perform the Wilcoxon test
    wilcox_result <- wilcox.test(formula, dat = dat)
    # Extract variable names from the formula
    response_var <- all.vars(formula)[1]
    group_var <- all.vars(formula)[2]
    # Extract medians for each group
    medians <- tapply(dat[[response_var]], dat[[group_var]], median, na.rm = TRUE)
    # Return a list with medians and the p-value
    return(list(medians = medians, p_value = wilcox_result$p.value))
}

# Helper function to generate and save violin plot
generate_violin_plot <- function(dat, feature, comp, filename, group, split = NULL, p_value = NULL, custom_colors, geneset_name = NULL) {
  violin_plot <- VlnPlot(dat, features = feature, pt.size = 0, group.by = group, split.by = split, cols = custom_colors)
  title_text <- paste0(geneset_name, "\n", comp, " (Wilcoxon p = ", format(p_value, scientific = TRUE), ")")
  violin_plot <- violin_plot + ggtitle(title_text) + theme(plot.title = element_text(size = 10))
  ggsave(filename = filename, plot = violin_plot, width = 7, height = 7)
}


# Main function
perform_wilcox_and_plots <- function(features, dat, name, geneset_names) {
# Initialize a data frame for storing results
  results_df <- data.frame(Feature=character(), Comparison=character(), Group1=character(),
                           Group2=character(), Median1=numeric(), Median2=numeric(), P_Value=numeric(), stringsAsFactors=FALSE)

  for (i in seq_along(features)) {
  feature <- features[i]
  geneset_name <- geneset_names[i]
    # Data preparation remains the same
    data_for_tests <- FetchData(dat, vars = c(feature, "MES", "sample_type"))
    colnames(data_for_tests) <- c("score", "MES", "sample_type")

    # MES vs non-MES
    comp <- "MES vs non-MES"
    mes_result <- perform_wilcox_test(data_for_tests, as.formula('score ~ MES'))
    generate_violin_plot(dat, feature, comp, paste0(feature, "_MES.png"), group = "MES", custom_colors = c("#387ee0", "gray"), p_value = mes_result$p_value, geneset_name = geneset_name)
    results_df <- rbind(results_df, data.frame(Feature=feature, Comparison=comp, Group1="MES", Group2="non-MES",
                                                Median1=mes_result$medians[1], Median2=mes_result$medians[2], P_Value=mes_result$p_value))
    
    # Primary vs Recurrent
    comp <- "Primary vs Recurrent"
    sample_type_result <- perform_wilcox_test(data_for_tests, as.formula('score ~ sample_type'))
    generate_violin_plot(dat, feature, comp, paste0(feature, "_sample_type.png"), group = "sample_type", custom_colors = c("#68F6B4", "#F58D9F"), p_value = sample_type_result$p_value, geneset_name = geneset_name)
    results_df <- rbind(results_df, data.frame(Feature=feature, Comparison=comp, Group1="Primary", Group2="Recurrent",
                                                Median1=sample_type_result$medians[1], Median2=sample_type_result$medians[2], P_Value=sample_type_result$p_value))
    
    # MES vs non-MES within each sample_type category
    for (sample_type in unique(dat$sample_type)) {
      comp <- paste0("MES vs non-MES within ", sample_type)
      dat_sample_type <- subset(dat, sample_type == sample_type)
      data_subset_for_tests <- data_for_tests[data_for_tests$sample_type == sample_type,]
      sample_type_mes_result <- perform_wilcox_test(data_subset_for_tests, as.formula('score ~ MES'))
      results_df <- rbind(results_df, data.frame(Feature=feature, Comparison=comp, Group1="MES", Group2="non-MES",
                                                  Median1=sample_type_mes_result$medians[1], Median2=sample_type_mes_result$medians[2], P_Value=sample_type_mes_result$p_value))
    }
    # plot split.by violing plot
    comp <- "MES vs non-MES within sample_type"
    generate_violin_plot(dat, feature, comp, paste0(feature, "_split_sample_type.png"), group = "MES", split = "sample_type", custom_colors = c("#68F6B4", "#F58D9F"), geneset_name = geneset_name)

    # Primary vs Recurrent within each MES category
      for (mes in unique(dat$MES)) {
        comp <- paste0("Primary vs Recurrent within ", mes)
        dat_mes <- subset(dat, MES == mes)
        data_subset_for_tests <- data_for_tests[data_for_tests$MES == mes,]
        mes_sample_type_result <- perform_wilcox_test(data_subset_for_tests, as.formula('score ~ sample_type'))
        results_df <- rbind(results_df, data.frame(Feature=feature, Comparison=comp, Group1="Primary", Group2="Recurrent",
                                                  Median1=mes_sample_type_result$medians[1], Median2=mes_sample_type_result$medians[2], P_Value=mes_sample_type_result$p_value))
      }
    # plot split.by violing plot
    comp <- "Primary vs Recurrent within MES"
    generate_violin_plot(dat, feature, comp, paste0(feature, "_split_MES.png"), group = "sample_type", split = "MES", custom_colors = c("#387ee0", "gray"), geneset_name = geneset_name)

  }
  # Save results to CSV
  write.csv(results_df, paste0(name, "_all_comparisons.csv"), row.names = FALSE)
}

# Define the 5 gene sets of interest
geneset_names <- c(
  "BIOCARTA_MAPK_PATHWAY",
  "BIOCARTA_ERK_PATHWAY",
  "KEGG_MAPK_SIGNALING_PATHWAY",
  "REACTOME_RAF_ACTIVATION",
  "REACTOME_ERK_MAPK_TARGETS",
  "REACTOME_MAPK1_ERK2_ACTIVATION"
)

# Pull msigdbr data for human
msigdb_data <- msigdbr(species = "Homo sapiens")

# Filter for selected gene sets
filtered_sets <- msigdb_data %>% filter(gs_name %in% geneset_names)

# Split into named vectors
geneset_list <- split(filtered_sets$gene_symbol, filtered_sets$gs_name)
geneset_list[["MAPK_MPAS"]] <- c(
  "PHLDA1", "SPRY2", "SPRY4", "DUSP4", "DUSP6",
  "CCND1", "EPHA2", "EPHA4", "ETV4", "ETV5"
)
mek_activation <- read_excel("ref/Dry-Pratilas-Signature.xlsx")
mek_activation <- as.vector(mek_activation$gene)
geneset_list[["MEK_ACTIVATION"]] <- mek_activation

for (name in names(geneset_list)) {
   # Add Module Score
  dat <- AddModuleScore(dat, features = list(geneset_list[[name]]), name = name)
  }

# Generate the feature names from geneset_list keys and call the function
feature_names <- paste0(names(geneset_list), "1")
perform_wilcox_and_plots(features = feature_names, dat = dat, name = "genesets", geneset_names = names(geneset_list))

##########################################
### ---------- genesets overlap ---------- ###
##########################################
# Create a named vector of gene set names
geneset_list[["MEK_ACTIVATION"]] <- mek_activation
# Load required packages
library(ComplexHeatmap)
library(circlize)

# List of gene set names
geneset_names_all <- names(geneset_list)

# Create empty matrix for Jaccard similarity
jaccard_matrix <- matrix(
  0, 
  nrow = length(geneset_names_all), 
  ncol = length(geneset_names_all),
  dimnames = list(geneset_names_all, geneset_names_all)
)

# Fill in Jaccard index values
for (i in geneset_names_all) {
  for (j in geneset_names_all) {
    set_i <- geneset_list[[i]]
    set_j <- geneset_list[[j]]
    intersection <- length(intersect(set_i, set_j))
    union <- length(union(set_i, set_j))
    jaccard_matrix[i, j] <- intersection / union
  }
}

# Plot heatmap using ComplexHeatmap
jaccard_heatmap<- Heatmap(jaccard_matrix,
        name = "Jaccard Index",
        col = colorRamp2(c(0, 0.5, 1), c("white", "skyblue", "darkblue")),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(title = "Similarity"))
# Save as PNG
png("genesets_jaccard_heatmap.png", width = 6, height = 5, units = "in", res = 300)
draw(jaccard_heatmap)
dev.off()

library(dplyr)
library(tibble)
library(purrr)

# Get all pairwise combinations of gene set names
geneset_combos <- expand.grid(Set1 = geneset_names_all, Set2 = geneset_names_all, stringsAsFactors = FALSE)

# Calculate overlap details along with individual set sizes
overlap_details <- geneset_combos %>%
  rowwise() %>%
  mutate(
    Set1_Size = length(geneset_list[[Set1]]),
    Set2_Size = length(geneset_list[[Set2]]),
    Overlapping_Genes = list(intersect(geneset_list[[Set1]], geneset_list[[Set2]])),
    Overlap_Count = length(Overlapping_Genes),
    Overlap_Genes = paste(Overlapping_Genes, collapse = ", ")
  ) %>%
  ungroup() %>%
  select(Set1, Set2, Set1_Size, Set2_Size, Overlap_Count, Overlap_Genes)

# View the resulting data frame
print(overlap_details)
# Save to CSV
write.csv(overlap_details, "geneset_overlap_details.csv", row.names = FALSE)

##################################
### ---------- dotplot  ---------- ###
##################################
mek_activation <- read_excel("ref/Dry-Pratilas-Signature.xlsx")
mek_activation <- as.vector(mek_activation$gene)

feature_dir <- "neftel/MES/dotplot"
ifelse(!dir.exists(file.path(main_dir, feature_dir)), 
       dir.create(file.path(main_dir, feature_dir)), FALSE)
setwd(file.path(main_dir, feature_dir))

#by MES, split by MES
Idents(dat)= "MES"

# Define gene sets to exclude
excluded_sets <- c("BIOCARTA_MAPK_PATHWAY", "KEGG_MAPK_SIGNALING_PATHWAY")

# Filter geneset list
geneset_list_filtered <- geneset_list[!names(geneset_list) %in% excluded_sets]

# Generate dotplots
for (set_name in names(geneset_list_filtered)) {
  
  feature_genes <- geneset_list_filtered[[set_name]]
  
  # Keep only genes present in data and remove duplicates
  feature_genes <- unique(feature_genes[feature_genes %in% rownames(dat)])
  
  if (length(feature_genes) == 0) next  # skip empty sets

  dotplot <- DotPlot(dat, features = feature_genes, cols = c("lightgray", "#FF0000")) +
    RotatedAxis() +
    guides(color = guide_colorbar(title = "Average Expression")) +
    ggtitle(paste("DotPlot:", set_name))
  
  ggsave(
    filename = paste0(set_name, "_dotplot_MES.png"),
    plot = dotplot,
    width = 10,
    height = 4)
}
