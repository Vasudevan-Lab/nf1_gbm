library(Seurat)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(ggpointdensity)
library(RColorBrewer)
library(stringr)
library(nichenetr)

#load data
dat <- readRDS("./rds/sb28_xeno_integrated_annotated_tumor_only.rds")
# Remove all columns ending with '_Score1' from meta.data
dat@meta.data <- dat@meta.data[, !grepl("_Score1$", colnames(dat@meta.data))]
# remove pano_ column
dat@meta.data <- dat@meta.data[, !grepl("pano_", colnames(dat@meta.data))]
# add new column called sample to rename orig.ident
dat@meta.data$sample <- dat@meta.data$orig.ident
tumor <- dat

# add metadata based on Neftel et al., 2019
reference_directory <- "ref/neftel_reference_sheets"
resultsdirectory <- "/sb28_10x/sel-ctrl-tumor/neftel"

#Calculate a MES, AC, OPC and NPC-like score for each cell in the dataset based on the expression of genes provided by Neftel
## Load the Neftel gene expression signatures
signatures <- read.table(file = paste0(reference_directory, "/neftel.IDHwt.GBM.MetaModules.tsv", sep = ""), header = TRUE)

##For MESlike and NPClike features, we extract the unique signatures from both MES/NPClike1 and 2. We then remove any NA values. 
MESlike_features <- unique(c(signatures$MESlike2, signatures$MESlike1)) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)]

NPClike_features <- unique(c(signatures$NPClike2, signatures$NPClike1)) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)]

# Adding G1.S and G2.M to the signatures list
signatures_list <- list(
  "MESlike" = as.character(MESlike_features),
  "AClike" = as.character(signatures$AClike) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)],
  "OPClike" = as.character(signatures$OPClike) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)],
  "NPClike" = as.character(NPClike_features),
  "MESlike1" = as.character(signatures$MESlike1) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)],
  "MESlike2" = as.character(signatures$MESlike2) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)],
  "NPClike1" = as.character(signatures$NPClike1) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)],
  "NPClike2" = as.character(signatures$NPClike2) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)],
  "G1.S" = as.character(signatures$G1.S) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)],
  "G2.M" = as.character(signatures$G2.M) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)]
)

signatures_list <- lapply(signatures_list, function(x){
  x[!is.na(x)]
})

## Calculate MESlike, AClike, OPClike and NPClike scores for each sample ####
tumor <- AddModuleScore(
  object = tumor,
  features = signatures_list[1:4],
  name = c("MESlike", "AClike", "OPClike", "NPClike"),
  search = TRUE
)

# Rename the columns in the metadata.
old_names <- c("MESlike1", "AClike2", "OPClike3", "NPClike4")
new_names <- c("MESlike.Score", "AClike.Score", "OPClike.Score", "NPClike.Score")
names(tumor@meta.data)[names(tumor@meta.data) %in% old_names] <- new_names

###We now run some calculations before generating the Neftel_Quadrants. 
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n] # a function that returns the position of n-th largest

## Calculate max scores for OPC-NPC vs AC-MES
##Here, the function assesses the  NPC and OPC-like score and takes the largest value of the two and creates a new column with the largest value.
tumor@meta.data$MaxScore_OPC.NPC <- apply(tumor@meta.data[, c("OPClike.Score", "NPClike.Score")], 1, function(x) x[which.max(x)])
tumor@meta.data$MaxScore_AC.MES <- apply(tumor@meta.data[, c("AClike.Score", "MESlike.Score")], 1, function(x) x[which.max(x)])
## Calculate D1. 
## RATIONALE: If for a CELL, MaxScore_OPC.NPC is the value of the OPC score and MaxScore_AC.MES is the value of the AC score, and IF OPC score > AC score then D1 (MaxScore_OPC.NPC - MaxScore_AC.MES) would be > 0. The CELL would be more OPClike than AClike. IF OPC score < AC score then D1 (MaxScore_OPC.NPC - MaxScore_AC.MES) would be < 0. The CELL would be more AClike than OPClike.
tumor@meta.data$D1 <- with(tumor@meta.data, MaxScore_OPC.NPC - MaxScore_AC.MES)
## Classify D1 > 0 as OPC-NPC cells, and D1 < 0 as AC-MES cells
tumor@meta.data$combGroup1 <- c("EMPTY")
tumor@meta.data$combGroup1[tumor@meta.data$D1 > 0] <- "OPC_NPC"
tumor@meta.data$combGroup1[tumor@meta.data$D1 < 0] <- "AC_MES"

## Calculate max scores for AC-OPC vs MES-NPC
##Here, the function assesses the  AC and OPC-like score and takes the largest value of the two and creates a new column with the largest value.
tumor@meta.data$MaxScore_AC.OPC <- apply(tumor@meta.data[, c("AClike.Score", "OPClike.Score")], 1, function(x) x[which.max(x)])
tumor@meta.data$MaxScore_MES.NPC <- apply(tumor@meta.data[, c("MESlike.Score", "NPClike.Score")], 1, function(x) x[which.max(x)])

## Calculate D2
## RATIONALE: If for a CELL, MaxScore_AC.OPC is the value of the OPC score and MaxScore_MES.NPC is the value of the MES score, and IF OPC score > MES score then D2 (MaxScore_AC.OPC - MaxScore_MES.NPC) would be > 0. The CELL would be more OPClike than MESlike. IF OPC score < AC score then D2 (MaxScore_AC.OPC - MaxScore_MES.NPC) would be < 0. The CELL would be more MESlike than OPClike.
tumor@meta.data$D2 <- with(tumor@meta.data, MaxScore_AC.OPC - MaxScore_MES.NPC)
## Classify D2 > 0 as AC-OPC cells, and D1 < 0 as MES-NPC cells
tumor@meta.data$combGroup2 <- c("EMPTY")
tumor@meta.data$combGroup2[tumor@meta.data$D2 > 0] <- "AC_OPC"
tumor@meta.data$combGroup2[tumor@meta.data$D2 < 0] <- "MES_NPC"

##We will now calculate a relative metamodule x and y score. 
##We want NPC and OPC-like cells to have a positive y-value and AC and MES-like cells to have a negative y-value.
##For cells that are classified as 'OPC_NPC' we log2 transform the difference between OPClike and NPClike scores.
##For cells that are classified as 'OPC_NPC' we (-)log2 transform the difference between OPClike and NPClike scores. A positive value will indicate OPC/NPC and a negative value AC/MES.
tumor@meta.data$rel.metamodule.y.score[tumor@meta.data$combGroup1 == "OPC_NPC"] <- with(tumor@meta.data[tumor@meta.data$combGroup1 == "OPC_NPC",], log2((abs(OPClike.Score - NPClike.Score)) + 1))
tumor@meta.data$rel.metamodule.y.score[tumor@meta.data$combGroup1 == "AC_MES"] <- with(tumor@meta.data[tumor@meta.data$combGroup1 == "AC_MES",], -log2((abs(AClike.Score - MESlike.Score)) + 1)) ## negative number for AC-MES

##We want NPC and MES-like cells to have a positive x-value and AC and OPC-like cells to have a negative y-value.
##For cells that are classified as 'AC_OPC' we -log2 transform the difference between OPClike and AClike scores.
##For cells that are classified as 'MES_NPC' we log2 transform the difference between NPClike and MESClike scores. A positive value will indicate MES/NPC and a negative value AC/OPC.
tumor@meta.data$rel.metamodule.x.score[tumor@meta.data$combGroup2 == "AC_OPC"] <- with(tumor@meta.data[tumor@meta.data$combGroup2 == "AC_OPC",], -log2((abs(AClike.Score - OPClike.Score)) + 1)) ## negative number for AC-OPC
tumor@meta.data$rel.metamodule.x.score[tumor@meta.data$combGroup2 == "MES_NPC"] <- with(tumor@meta.data[tumor@meta.data$combGroup2 == "MES_NPC",], log2((abs(MESlike.Score - NPClike.Score)) + 1))

## We can now classify cells into different cell states.
## Left-Top Quadrant: OPC-like - Negative x-value and Positive y-value.
## Right-Top Quadrant: NPC-like - Positive x-value and Positive y-value.
## Left-Bottom Quadrant: AC-like - Negative x-value and Negative y-value.
## Right-Bottom Quadrant: MES-like - Positive x-value and Negative y-value.

## Add new column in metadata with cell-state name.
tumor@meta.data$Quadrant <- c("EMPTY")
tumor@meta.data$Quadrant[tumor@meta.data$rel.metamodule.x.score < 0 & tumor@meta.data$rel.metamodule.y.score < 0] <- "AClike"
tumor@meta.data$Quadrant[tumor@meta.data$rel.metamodule.x.score < 0 & tumor@meta.data$rel.metamodule.y.score > 0] <- "OPClike"
tumor@meta.data$Quadrant[tumor@meta.data$rel.metamodule.x.score > 0 & tumor@meta.data$rel.metamodule.y.score > 0] <- "NPClike"
tumor@meta.data$Quadrant[tumor@meta.data$rel.metamodule.x.score > 0 & tumor@meta.data$rel.metamodule.y.score < 0] <- "MESlike"

##Calculate the number of cells in each quadrant; per patient
cellstate_perpatient <- tumor@meta.data %>% group_by(sample, Quadrant) %>% tally()
write.csv(as.data.frame(cellstate_perpatient), file = paste0(resultsdirectory, "/GBM_All_sheets/GBM_cellstate_perpatient.csv"), quote = FALSE)

####We now Calculate MESlike1,2 and NPClike1,2 scores for each sample #####
tumor <- AddModuleScore(
  object = tumor,
  features = signatures_list[5:8],
  name = c("MESlike1", "MESlike2", "NPClike1", "NPClike2"),
  search = TRUE
)

##Rename the columns in the metadata.
old_names <- c("MESlike11", "MESlike22", "NPClike13", "NPClike24")
new_names <- c("MESlike1.Score", "MESlike2.Score", "NPClike1.Score", "NPClike2.Score")
names(tumor@meta.data)[names(tumor@meta.data) %in% old_names] <- new_names

##We create new metadata columns to classify cells as MESlike1, 2 or not MESlike. 
##We do this ONLY for cells that we classified as MESlike above. 
tumor$MESlike <- ifelse(tumor$MESlike1.Score > tumor$MESlike2.Score & tumor$Quadrant == "MESlike", "MESlike1", ifelse(tumor$MESlike1.Score < tumor$MESlike2.Score & tumor$Quadrant == "MESlike", "MESlike2", "notMESlike"))
unique(tumor$MESlike)

##We create new metadata columns to classify cells as NPClike1, 2 or not NPClike. ####
##We do this ONLY for cells that we classified as NPClike above. 
tumor$NPClike <- ifelse(tumor$NPClike1.Score > tumor$NPClike2.Score & tumor$Quadrant == "NPClike", "NPClike1", 
       ifelse(tumor$NPClike1.Score < tumor$NPClike2.Score & tumor$Quadrant == "NPClike", "NPClike2", "notNPClike"))
unique(tumor$NPClike)

# We will also calculate cycling scores.

tumor <- AddModuleScore(
  object = tumor,
  features = signatures_list[9:10],
  name = c("G1.S", "G2.M"),
  search = TRUE
)

##Rename the columns in the metadata.
old_names <- c("G1.S1", "G2.M2")
new_names <- c("G1.S.Score", "G2.M.Score")
names(tumor@meta.data)[names(tumor@meta.data) %in% old_names] <- new_names

# Fit a normal distribution to the G1S and G2M scores in the 'tumor' dataset
fit_g1s <- fitdistr(tumor$G1.S.Score, "normal")
fit_g2m <- fitdistr(tumor$G2.M.Score, "normal")

# Define a threshold for p < 0.001 based on the fitted normal distribution
threshold_g1s <- qnorm(0.999, mean = fit_g1s$estimate["mean"], sd = fit_g1s$estimate["sd"])
threshold_g2m <- qnorm(0.999, mean = fit_g2m$estimate["mean"], sd = fit_g2m$estimate["sd"])

# Classify cells as cycling or non-cycling based on the thresholds
tumor$Cell_Cycle <- "NonCycling"
tumor$Cell_Cycle[tumor$G1.S.Score > threshold_g1s | tumor$G2.M.Score > threshold_g2m] <- "Cycling"
unique(tumor$Cell_Cycle)

###Count cycling and non-cycling cells; per quadrant, per condition####
cellstate_cellcycle <- tumor@meta.data %>% group_by(Cell_Cycle, Quadrant) %>% tally()
write.csv(as.data.frame(cellstate_cellcycle), file = paste0(resultsdirectory, "/GBM_All_sheets/GBM_cellstate_cellcycle.csv"), quote = FALSE)

##Save MetaData and Seurat Object.
write.csv(tumor@meta.data, file = paste0(resultsdirectory,"/tumor_MetaData.csv"))
# save(tumor, file = paste0(resultsdirectory,"/tumor_All.Rdata"))

# plot UMAP based on Quadrant
TumorColors <- c(
  "MESlike" = "#3385BB",  # blue
  "AClike" = "#6DBD57",   # green
  "OPClike" = "#FE8D19",  # orange
  "NPClike" = "#E42622")   # red
plot_umap <- DimPlot(tumor, reduction = "umap", group.by = "Quadrant", cols = TumorColors)
ggsave(plot_umap, filename = paste0(resultsdirectory,"/plots/GBM_umap_quadrant.pdf"), height = 5, width = 6, dpi = 300, useDingbats = FALSE)


## FIGURE: Recreate Neftel meta-module plots with ggplot2. ####

# 1) plot by sample

# Calculate the percentage of cells per quadrant and per sample
percentage_data <- tumor@meta.data %>%
  group_by(sample, Quadrant) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100)

percentage_data <- percentage_data %>%
  ungroup() %>%
  complete(sample, Quadrant = c("AClike", "MESlike", "NPClike", "OPClike"), fill = list(count = 0, percentage = 0))


# Create label_data with type, label, and correct x/y coordinates for positioning
label_data <- percentage_data %>%
  mutate(type = str_replace(Quadrant, "like", ""),  # Remove 'like' suffix to get the clean type
         x = case_when(
           type == "OPC" ~ -1,  # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ 1    # MES -> Lower right
         ),
         y = case_when(
           type == "OPC" ~ 1,   # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ -1   # MES -> Lower right
         ),
         label = paste0(type, ": ", round(percentage, 1), "%"))

int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0] ## function to generate only whole-integer breaks for x and y axis ticks
scatterPlot <- ggplot(tumor@meta.data, aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score)) + ## transform to change the order of the plots in facet_wrap()
  geom_vline(xintercept=0, linetype="dotted", color = "black", size = 0.5) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size = 0.5) +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_gradientn(
    colours = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(25),
    limits = c(0, 2000),  # Set limits to 0-2000
    breaks = c(0, 2000),  # Only show 0 and 2000 in the legend
    labels = c("0", "2000")) +
  guides(color = guide_colorbar(
    title = "Density",  # Add "Density" as the title of the legend
    title.position = "top",  # Position the title at the top
    title.hjust = 0.5,
    override.aes = list(size = 5),
    barwidth = 10,
    barheight = 1)) +
  facet_wrap( ~ sample) + ##Condition or Condition_PatientID to generate split scatter plots by Condition or PatientID
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),         
    legend.title = element_text(size = 20),       
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20),
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5)
  ) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) + ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25))+
  geom_text(data = label_data, aes(x = x, y = y, label = label),
            size = 5, inherit.aes = FALSE, hjust = ifelse(label_data$x == 1, 1, 0), vjust = ifelse(label_data$y == 1, 1, 0))

ggsave(scatterPlot, filename = paste0(resultsdirectory,"/plots/GBM_metamodules_bysample_scale0-2000.pdf"), height = 8, width = 12, dpi = 300, useDingbats = FALSE)

# 2) plot by stim
# Calculate the percentage of cells per quadrant and per stim
percentage_data <- tumor@meta.data %>%
  group_by(stim, Quadrant) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100)

percentage_data <- percentage_data %>%
  ungroup() %>%
  complete(stim, Quadrant = c("AClike", "MESlike", "NPClike", "OPClike"), fill = list(count = 0, percentage = 0))


# Create label_data with type, label, and correct x/y coordinates for positioning
label_data <- percentage_data %>%
  mutate(type = str_replace(Quadrant, "like", ""),  # Remove 'like' suffix to get the clean type
         x = case_when(
           type == "OPC" ~ -1,  # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ 1    # MES -> Lower right
         ),
         y = case_when(
           type == "OPC" ~ 1,   # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ -1   # MES -> Lower right
         ),
         label = paste0(type, ": ", round(percentage, 1), "%"))

int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0] ## function to generate only whole-integer breaks for x and y axis ticks
scatterPlot <- ggplot(tumor@meta.data, aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score)) + ## transform to change the order of the plots in facet_wrap()
  geom_vline(xintercept=0, linetype="dotted", color = "black", size = 0.5) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size = 0.5) +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_gradientn(
    colours = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(25),
    limits = c(0, 2000),  # Set limits to 0-2000
    breaks = c(0, 2000),  # Only show 0 and 2000 in the legend
    labels = c("0", "2000")) +
  guides(color = guide_colorbar(
    title = "Density",  # Add "Density" as the title of the legend
    title.position = "top",  # Position the title at the top
    title.hjust = 0.5,
    barwidth = 10,
    barheight = 1)) +
  facet_wrap( ~ stim) + ##Condition or Condition_PatientID to generate split scatter plots by Condition or PatientID
  theme_bw() +
  theme(legend.position = "top", axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), strip.text = element_text(size = 20), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) + ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25))+
  geom_text(data = label_data, aes(x = x, y = y, label = label),
            size = 5, inherit.aes = FALSE, hjust = ifelse(label_data$x == 1, 1, 0), vjust = ifelse(label_data$y == 1, 1, 0))

ggsave(scatterPlot, filename = paste0(resultsdirectory,"/plots/GBM_metamodules_bystim_scale0-2000.pdf"), height = 8, width = 12, dpi = 300, useDingbats = FALSE)


# 3) plot combined

# Function to generate only whole-integer breaks for x and y axis ticks
int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]

# Calculate the percentage of cells per quadrant across all samples
percentage_data <- tumor@meta.data %>%
  group_by(Quadrant) %>%
  summarise(count = n()) %>%
  mutate(total = sum(count),
         percentage = (count / total) * 100) %>%
  ungroup() %>%
  complete(Quadrant = c("AClike", "MESlike", "NPClike", "OPClike"), fill = list(count = 0, percentage = 0))

# Create label_data with type, label, and correct x/y coordinates for positioning
label_data <- percentage_data %>%
  mutate(type = str_replace(Quadrant, "like", ""),  # Remove 'like' suffix to get the clean type
         x = case_when(
           type == "OPC" ~ -1,  # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ 1    # MES -> Lower right
         ),
         y = case_when(
           type == "OPC" ~ 1,   # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ -1   # MES -> Lower right
         ),
         label = paste0(type, ": ", round(percentage, 1), "%"))

# Create the combined scatter plot without faceting
scatterPlot <- ggplot(tumor@meta.data, aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score)) +
  geom_vline(xintercept=0, linetype="dotted", color = "black", size = 0.5) +
  geom_hline(yintercept=0, linetype="dotted", color = "black", size = 0.5) +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_gradientn(
    colours = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(25),
    limits = c(0, 5000),  # Set limits to 0-5000
    breaks = c(0, 5000),  # Only show 0 and 5000 in the legend
    labels = c("0", "5000")) +
  guides(color = guide_colorbar(
    title = "Density",  # Add "Density" as the title of the legend
    title.position = "top",  # Position the title at the top
    title.hjust = 0.5,
    barwidth = 10,
    barheight = 1)) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  geom_text(data = label_data, aes(x = x, y = y, label = label),
            size = 5, inherit.aes = FALSE,
            hjust = ifelse(label_data$x == 1, 1, 0),
            vjust = ifelse(label_data$y == 1, 1, 0))

# Save the plot to a file
ggsave(scatterPlot, filename = paste0(resultsdirectory,"/plots/GBM_metamodules_combined_scale0-5000.pdf"), height = 6, width = 5, dpi = 50, useDingbats = FALSE)

# 3) plot by stim for cycling cells
# Updated code for plotting cycling percentage per quadrant and stim condition

# Function to generate only whole-integer breaks for x and y axis ticks
int_breaks <- function(x, n = 5) pretty(x, n)[pretty(x, n) %% 1 == 0]

# Calculate cycling percentage per quadrant and stim, as before
cycling_percentage_data <- tumor@meta.data %>%
  group_by(stim, Quadrant) %>%
  summarise(total_count = n(),  # Total cell count per quadrant and stim
            cycling_count = sum(Cell_Cycle == "Cycling"),  # Cycling cell count per quadrant and stim
            .groups = 'drop') %>%
  mutate(percentage = (cycling_count / total_count) * 100) %>%
  complete(stim, Quadrant = c("AClike", "MESlike", "NPClike", "OPClike"), 
           fill = list(total_count = 0, cycling_count = 0, percentage = 0))

# Generate label data for cycling percentage in each quadrant
cycling_label_data <- cycling_percentage_data %>%
  mutate(type = str_replace(Quadrant, "like", ""),
         x = case_when(
           type == "OPC" ~ -1,  # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ 1    # MES -> Lower right
         ),
         y = case_when(
           type == "OPC" ~ 1,   # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ -1   # MES -> Lower right
         ),
         label = paste0(type, ": ", round(percentage, 1), "%"))

# Create scatter plot with cycling percentage labels instead of density
scatterPlot_cycling_percentage_by_stim <- ggplot(tumor@meta.data, aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  # Non-cycling cells as individual points for visual context
  geom_point(data = subset(tumor@meta.data, Cell_Cycle == "NonCycling"), color = "grey", size = 1.5, alpha = 0.5) +
  # Cycling cells as points (adjust color if needed for better contrast)
  geom_point(data = subset(tumor@meta.data, Cell_Cycle == "Cycling"), color = "orange", size = 1.5, alpha = 0.7) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5)) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  facet_wrap(~ stim) +  # Split scatter plots by stim condition
  # Add cycling percentage labels at predefined coordinates per quadrant
  geom_text(data = cycling_label_data, aes(x = x, y = y, label = label),
            size = 5, inherit.aes = FALSE,
            hjust = ifelse(cycling_label_data$x == 1, 1, 0),
            vjust = ifelse(cycling_label_data$y == 1, 1, 0))

# Save the plot
ggsave(scatterPlot_cycling_percentage_by_stim, filename = paste0(resultsdirectory, "/plots/GBM_metamodules_bystim_cycling_percentage.pdf"), height = 6, width = 10, dpi = 300, useDingbats = FALSE)

# 3) plot cycling based on G2.M

# Calculate cycling percentage per quadrant, excluding stim grouping
cycling_percentage_data <- tumor@meta.data %>%
  group_by(Quadrant) %>%
  summarise(total_count = n(),  # Total cell count per quadrant
            cycling_count = sum(Cell_Cycle == "Cycling"),  # Cycling cell count per quadrant
            .groups = 'drop') %>%
  mutate(percentage = (cycling_count / total_count) * 100) %>%
  complete(Quadrant = c("AClike", "MESlike", "NPClike", "OPClike"), 
           fill = list(total_count = 0, cycling_count = 0, percentage = 0))

# Generate label data for cycling percentage in each quadrant
cycling_label_data <- cycling_percentage_data %>%
  mutate(type = str_replace(Quadrant, "like", ""),
         x = case_when(
           type == "OPC" ~ -1,  # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ 1    # MES -> Lower right
         ),
         y = case_when(
           type == "OPC" ~ 1,   # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ -1   # MES -> Lower right
         ),
         label = paste0(type, ": ", round(percentage, 1), "%"))

# Create combined scatter plot with cycling percentage labels instead of density
scatterPlot_combined_cycling_percentage <- ggplot(tumor@meta.data, aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  # Non-cycling cells as individual points for visual context
  geom_point(data = subset(tumor@meta.data, Cell_Cycle == "NonCycling"), color = "grey", size = 1.5, alpha = 0.5) +
  # Cycling cells as points (adjust color if needed for better contrast)
  geom_point(data = subset(tumor@meta.data, Cell_Cycle == "Cycling"), color = "orange", size = 1.5, alpha = 0.7) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5)) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  # Add cycling percentage labels at predefined coordinates per quadrant
  geom_text(data = cycling_label_data, aes(x = x, y = y, label = label),
            size = 5, inherit.aes = FALSE,
            hjust = ifelse(cycling_label_data$x == 1, 1, 0),
            vjust = ifelse(cycling_label_data$y == 1, 1, 0))

# Save the combined plot
ggsave(scatterPlot_combined_cycling_percentage, filename = paste0(resultsdirectory, "/plots/GBM_metamodules_combined_cycling_percentage.pdf"), height = 5, width = 5, dpi = 50, useDingbats = FALSE)


# plot cycling based on G1.S and G2.M
# Calculate the percentage of cycling cells
cycling_percentage <- (sum(tumor@meta.data$Cell_Cycle == "Cycling") / nrow(tumor@meta.data)) * 100

# Create a scatter plot for G1/S vs G2/M scores
cell_cycle_scatter <- ggplot(tumor@meta.data, aes(x = G1.S.Score, y = G2.M.Score)) +
  geom_point(aes(color = Cell_Cycle), size = 2) +  # Plot points colored by cell cycle phase
  scale_color_manual(values = c("NonCycling" = "gray", "Cycling" = "orange")) +  # Assign colors for cycling and non-cycling
  # Add cycling percentage label at the upper-right corner
  annotate("text", x = Inf, y = Inf, label = paste0("Cycling: ", round(cycling_percentage, 1), "%"),
           hjust = 1.2, vjust = 2, size = 5, color = "black") +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5)) +
  xlab("G1/S score") +
  ylab("G2/M score") +
  guides(color = guide_legend(title = "Cell Cycle Phase"))

# Save the plot to a file
ggsave(cell_cycle_scatter, filename = paste0(resultsdirectory, "/plots/GBM_cell_cycle_G1S_G2M.pdf"), height = 5, width = 5, dpi = 50, useDingbats = FALSE)

# Calculate the percentage of cycling cells by stim
cycling_percentage_by_stim <- tumor@meta.data %>%
  group_by(stim) %>%
  summarise(cycling_percentage = (sum(Cell_Cycle == "Cycling") / n()) * 100)

# Create a scatter plot for G1/S vs G2/M scores, split by stim with cycling percentage labels
cell_cycle_scatter <- ggplot(tumor@meta.data, aes(x = G1.S.Score, y = G2.M.Score)) +
  geom_point(aes(color = Cell_Cycle), size = 2) +  # Plot points colored by cell cycle phase
  scale_color_manual(values = c("NonCycling" = "gray", "Cycling" = "orange")) +  # Assign colors for cycling and non-cycling
  facet_wrap(~ stim) +  # Split by stim condition
  # Add cycling percentage label for each stim at the upper-right corner of each facet
  geom_text(data = cycling_percentage_by_stim, 
            aes(x = Inf, y = Inf, label = paste0("Cycling: ", round(cycling_percentage, 1), "%")),
            hjust = 1.2, vjust = 2, size = 5, color = "black", inherit.aes = FALSE) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5)) +
  xlab("G1/S score") +
  ylab("G2/M score") +
  guides(color = guide_legend(title = "Cell Cycle Phase"))

# Save the plot to a file
ggsave(cell_cycle_scatter, filename = paste0(resultsdirectory, "/plots/GBM_cell_cycle_G1S_G2M_by_stim.pdf"), height = 6, width = 10, dpi = 300, useDingbats = FALSE)

# 5) plot by stim

# Extract the color corresponding to the density value of 0 from the color palette
color_palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(25)
color_zero_density <- color_palette[1]  # This is the color for the lowest density (0)

# Create a scatter plot showing density of cycling cells
scatterPlot_cycling_density <- ggplot(tumor@meta.data, aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.5) +  
  # Plot non-cycling cells as individual points, using the color for density 0
  geom_point(data = subset(tumor@meta.data, Cell_Cycle == "NonCycling"), color = color_zero_density, size = 2) +
  # Show density for cycling cells
  geom_pointdensity(data = subset(tumor@meta.data, Cell_Cycle == "Cycling"), adjust = 0.5) +
  scale_colour_gradientn(
    colours = color_palette,
    limits = c(0, 100),  # Adjust limits as needed based on your data
    breaks = c(0, 100),  # Show density range in the legend
    labels = c("0", "100")) +
  guides(color = guide_colorbar(
    title = "Cycling Cell Density",  # Legend title specific for cycling cells
    title.position = "top",
    title.hjust = 0.5,
    barwidth = 10,
    barheight = 1)) +
  facet_wrap( ~ stim) + ##Condition or Condition_PatientID to generate split scatter plots by Condition or PatientID
  theme_bw() +
  theme(legend.position = "top", axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), strip.text = element_text(size = 15), panel.border = element_rect(fill=NA, colour = "black", size=1.5)) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) + ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25))+
  geom_text(data = label_data, aes(x = x, y = y, label = label),
            size = 5, inherit.aes = FALSE, hjust = ifelse(label_data$x == 1, 1, 0), vjust = ifelse(label_data$y == 1, 1, 0))

# Save the plot to a file
ggsave(scatterPlot_cycling_density, filename = paste0(resultsdirectory, "/plots/GBM_metamodules_bystim_cycling_scale0-100.pdf"), height = 8, width = 12, dpi = 300, useDingbats = FALSE)
 
# plot cycling based on G1.S and G2.M
# Calculate cycling percentage for each stim condition
cycling_percentage_stim <- tumor@meta.data %>%
  group_by(stim) %>%
  summarize(cycling_percentage = (sum(Cell_Cycle == "Cycling") / n()) * 100)

# Create a scatter plot for G1/S vs G2/M scores, split by stim condition
cell_cycle_scatter <- ggplot(tumor@meta.data, aes(x = G1.S.Score, y = G2.M.Score)) +
  geom_point(aes(color = Cell_Cycle), size = 2) +  # Plot points colored by cell cycle phase
  scale_color_manual(values = c("NonCycling" = "gray", "Cycling" = "orange")) +  # Assign colors for cycling and non-cycling
  # Add cycling percentage label for each stim condition
  geom_text(data = cycling_percentage_stim,
            aes(x = Inf, y = Inf, label = paste0("Cycling: ", round(cycling_percentage, 1), "%")),
            hjust = 1.2, vjust = 2, size = 5, color = "black", inherit.aes = FALSE) +
  facet_wrap(~ stim) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5)) +
  xlab("G1/S score") +
  ylab("G2/M score") +
  guides(color = guide_legend(title = "Cell Cycle Phase"))
# Save the plot to a file
ggsave(cell_cycle_scatter, filename = paste0(resultsdirectory, "/plots/GBM_cell_cycle_G1S_G2M_bystim.pdf"), height = 5, width = 10, dpi = 300, useDingbats = FALSE)

# 4) plot by stim for cycling cells
# Updated code for plotting cycling percentage per quadrant and stim condition

# Calculate cycling percentage per quadrant and stim, as before
cycling_percentage_data <- tumor@meta.data %>%
  group_by(stim, Quadrant) %>%
  summarise(total_count = n(),  # Total cell count per quadrant and stim
            cycling_count = sum(Cell_Cycle == "Cycling"),  # Cycling cell count per quadrant and stim
            .groups = 'drop') %>%
  mutate(percentage = (cycling_count / total_count) * 100) %>%
  complete(stim, Quadrant = c("AClike", "MESlike", "NPClike", "OPClike"), 
           fill = list(total_count = 0, cycling_count = 0, percentage = 0))

# Generate label data for cycling percentage in each quadrant
cycling_label_data <- cycling_percentage_data %>%
  mutate(type = str_replace(Quadrant, "like", ""),
         x = case_when(
           type == "OPC" ~ -1,  # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ 1    # MES -> Lower right
         ),
         y = case_when(
           type == "OPC" ~ 1,   # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ -1   # MES -> Lower right
         ),
         label = paste0(type, ": ", round(percentage, 1), "%"))

# Create scatter plot with cycling percentage labels instead of density
scatterPlot_cycling_percentage_by_stim <- ggplot(tumor@meta.data, aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  # Non-cycling cells as individual points for visual context
  geom_point(data = subset(tumor@meta.data, Cell_Cycle == "NonCycling"), color = "grey", size = 1.5, alpha = 0.5) +
  # Cycling cells as points (adjust color if needed for better contrast)
  geom_point(data = subset(tumor@meta.data, Cell_Cycle == "Cycling"), color = "orange", size = 1.5, alpha = 0.7) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5)) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  facet_wrap(~ stim) +  # Split scatter plots by stim condition
  # Add cycling percentage labels at predefined coordinates per quadrant
  geom_text(data = cycling_label_data, aes(x = x, y = y, label = label),
            size = 5, inherit.aes = FALSE,
            hjust = ifelse(cycling_label_data$x == 1, 1, 0),
            vjust = ifelse(cycling_label_data$y == 1, 1, 0))

# Save the plot
ggsave(scatterPlot_cycling_percentage_by_stim, filename = paste0(resultsdirectory, "/plots/GBM_metamodules_bystim_cycling_percentage.pdf"), height = 6, width = 10, dpi = 300, useDingbats = FALSE)

cycling_percentage_data <- tumor@meta.data %>%
  group_by(sample, Quadrant) %>%
  summarise(total_count = n(),  # Total cell count per quadrant and stim
            cycling_count = sum(Cell_Cycle == "Cycling"),  # Cycling cell count per quadrant and stim
            .groups = 'drop') %>%
  mutate(percentage = (cycling_count / total_count) * 100) %>%
  complete(sample, Quadrant = c("AClike", "MESlike", "NPClike", "OPClike"), 
           fill = list(total_count = 0, cycling_count = 0, percentage = 0))

# Create scatter plot with cycling percentage labels instead of density
scatterPlot_cycling_percentage_by_sample <- ggplot(tumor@meta.data, aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  # Non-cycling cells as individual points for visual context
  geom_point(data = subset(tumor@meta.data, Cell_Cycle == "NonCycling"), color = "grey", size = 1.5, alpha = 0.5) +
  # Cycling cells as points (adjust color if needed for better contrast)
  geom_point(data = subset(tumor@meta.data, Cell_Cycle == "Cycling"), color = "orange", size = 1.5, alpha = 0.7) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        panel.border = element_rect(fill = NA, colour = "black", size = 1.5)) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  facet_wrap(~ sample) +  # Split scatter plots by stim condition
  # Add cycling percentage labels at predefined coordinates per quadrant
  geom_text(data = cycling_label_data, aes(x = x, y = y, label = label),
            size = 5, inherit.aes = FALSE,
            hjust = ifelse(cycling_label_data$x == 1, 1, 0),
            vjust = ifelse(cycling_label_data$y == 1, 1, 0))

# Save the plot
ggsave(scatterPlot_cycling_percentage_by_sample, filename = paste0(resultsdirectory, "/plots/GBM_metamodules_bysample_cycling_percentage.pdf"), height = 7, width = 10, dpi = 50, useDingbats = FALSE)



# 6) plot by sample, color by subclones (no NA)

# Load the subclone data
subclone_data <- readRDS("../sel-ctrl/SCEVAN_by_sample/combined_results.rds")

# Add subclone metadata to the tumor Seurat object
subclone_metadata <- subclone_data@meta.data[, "subclone", drop = FALSE]
cell_barcodes <- rownames(subclone_metadata)

# Ensure that only matching barcodes are used
matching_barcodes <- cell_barcodes[cell_barcodes %in% rownames(tumor@meta.data)]
tumor <- AddMetaData(object = tumor, metadata = subclone_metadata[matching_barcodes, , drop = FALSE])

# Create label_data with type, label, and correct x/y coordinates for positioning
label_data <- percentage_data %>%
  mutate(type = str_replace(Quadrant, "like", ""),  # Remove 'like' suffix to get the clean type
         x = case_when(
           type == "OPC" ~ -1,  # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ 1    # MES -> Lower right
         ),
         y = case_when(
           type == "OPC" ~ 1,   # OPC -> Upper left
           type == "NPC" ~ 1,   # NPC -> Upper right
           type == "AC" ~ -1,   # AC -> Lower left
           type == "MES" ~ -1   # MES -> Lower right
         ),
         label = type)


# Create scatter plot colored by subclones
scatterPlot_subclone <- ggplot(tumor@meta.data, aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score, color = factor(subclone))) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_manual(values = brewer.pal(8, 'Set1')) +
  guides(color = guide_legend(title = "Subclone", title.position = "top", title.hjust = 0.5), override.aes = list(size = 5)) +
  facet_wrap(~ sample) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),         
    legend.title = element_text(size = 20),       
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20),
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5)
  ) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  geom_text(data = label_data, aes(x = x, y = y, label = label), size = 5, inherit.aes = FALSE, hjust = 0.5, vjust = 0.5)

# Save the plot
ggsave(scatterPlot_subclone, filename = paste0(resultsdirectory, "/plots/GBM_metamodules_bysample_subclones.pdf"), height = 10, width = 12, dpi = 300, useDingbats = FALSE)


# Create scatter plot colored by subclones (no NA)
scatterPlot_subclone <- ggplot(tumor@meta.data[!is.na(tumor@meta.data$subclone), ], aes(x = rel.metamodule.x.score, y = rel.metamodule.y.score, color = factor(subclone))) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 0.5) +
  geom_point(size = 2) +
  geom_pointdensity(adjust = 0.5) +
  scale_colour_manual(values = brewer.pal(8, 'Set1')) +
  guides(color = guide_legend(title = "Subclone", title.position = "top", title.hjust = 0.5, override.aes = list(size = 5))) +
  facet_wrap(~ sample) +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 20),         
    legend.title = element_text(size = 20),       
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    strip.text = element_text(size = 20),
    panel.border = element_rect(fill = NA, colour = "black", size = 1.5)
  ) +
  xlab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  ylab(expression(atop("Relative meta-module score", paste("[log2(|SC1-SC2|+1)]")))) +
  scale_x_continuous(breaks = int_breaks) +
  scale_y_continuous(breaks = int_breaks) +
  coord_cartesian(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  geom_text(data = label_data, aes(x = x, y = y, label = label), size = 5, inherit.aes = FALSE, hjust = 0.5, vjust = 0.5)

# Save the plot
ggsave(scatterPlot_subclone, filename = paste0(resultsdirectory, "/plots/GBM_metamodules_bysample_subclones_noNA.pdf"), height = 10, width = 12, dpi = 300, useDingbats = FALSE)

library(dplyr)

# Calculate count of NA and non-NA subclone values per sample and add percentage columns
cell_count_summary <- tumor@meta.data %>%
  group_by(sample) %>%
  summarise(
    NA_count = sum(is.na(subclone)),
    non_NA_count = sum(!is.na(subclone)),
    total_count = n()
  ) %>%
  mutate(
    NA_percentage = (NA_count / total_count) * 100,
    non_NA_percentage = (non_NA_count / total_count) * 100
  )

# View the summary
print(cell_count_summary)

#####################################################
### --------------- update c-f --------------- ###
#####################################################
# 5c: plot UMAP based on Cell_Cycle
Cycling_colors <- c("Cycling" = "orange", "NonCycling" = "gray")
plot_umap <- DimPlot(tumor, reduction = "umap", group.by = "Cell_Cycle", cols = Cycling_colors)
ggsave(plot_umap, filename = paste0(resultsdirectory,"/plots/GBM_umap_cycling.pdf"), height = 5, width = 6, dpi = 300, useDingbats = FALSE)

# 5e: violin plot of CRISPRi genes
# Create a new column called 'MES' in the metadata of 'tumor'
tumor$MES <- ifelse(tumor$MESlike %in% c("MESlike1", "MESlike2"), "MES", "nonMES")


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
generate_violin_plot <- function(dat, feature, comp, filename, group, split = NULL, p_value = NULL, custom_colors) {
  violin_plot <- VlnPlot(dat, features = feature, pt.size = 0, group.by = group, split.by = split, cols = custom_colors)
  # add p value to title
  violin_plot <- violin_plot + ggtitle(paste0(comp, " Wilcox-test p-value: ", format(p_value, scientific = TRUE))) + theme(plot.title = element_text(size=10))
  # Save the plot
  ggsave(filename = filename, plot = violin_plot, width = 7, height = 7)
}

# Main function
perform_wilcox_and_plots <- function(features, dat, name) {
# Initialize a data frame for storing results
  results_df <- data.frame(Feature=character(), Comparison=character(), Group1=character(),
                           Group2=character(), Median1=numeric(), Median2=numeric(), P_Value=numeric(), stringsAsFactors=FALSE)

  for (feature in features) {
    # Data preparation remains the same
    data_for_tests <- FetchData(dat, vars = c(feature, "MES", "stim"))
    colnames(data_for_tests) <- c("score", "MES", "stim")

    # MES vs non-MES
    comp <- "MES vs non-MES"
    mes_result <- perform_wilcox_test(data_for_tests, as.formula('score ~ MES'))
    generate_violin_plot(dat, feature, comp, paste0(feature, "_MES.png"), group = "MES", custom_colors = c("#387ee0", "gray"), p_value = mes_result$p_value)
    results_df <- rbind(results_df, data.frame(Feature=feature, Comparison=comp, Group1="MES", Group2="non-MES",
                                                Median1=mes_result$medians[1], Median2=mes_result$medians[2], P_Value=mes_result$p_value))
    
    # Selumetinib vs Vehicle
    comp <- "Selumetinib vs Vehicle"
    stim_result <- perform_wilcox_test(data_for_tests, as.formula('score ~ stim'))
    generate_violin_plot(dat, feature, comp, paste0(feature, "_Stim.png"), group = "stim", custom_colors = c("#FF0000", "#0000FF"), p_value = stim_result$p_value)
    results_df <- rbind(results_df, data.frame(Feature=feature, Comparison=comp, Group1="Selumetinib", Group2="Vehicle",
                                                Median1=stim_result$medians[1], Median2=stim_result$medians[2], P_Value=stim_result$p_value))
    # MES vs non-MES within each stim category
    for (stim in unique(dat$stim)) {
      comp <- paste0("MES vs non-MES within ", stim)
      dat_stim <- subset(dat, stim == stim)
      data_subset_for_tests <- data_for_tests[data_for_tests$stim == stim,]
      stim_mes_result <- perform_wilcox_test(data_subset_for_tests, as.formula('score ~ MES'))
        generate_violin_plot(dat_stim, feature, comp, paste0(feature, "_", stim, ".png"), group = "MES", p_value = stim_mes_result$p_value, custom_colors = c("#387ee0", "gray"))
      results_df <- rbind(results_df, data.frame(Feature=feature, Comparison=comp, Group1="MES", Group2="non-MES",
                                                  Median1=stim_mes_result$medians[1], Median2=stim_mes_result$medians[2], P_Value=stim_mes_result$p_value))
    }
    # Selumetinib vs Vehicle within each MES category
      for (mes in unique(dat$MES)) {
        comp <- paste0("Selumetinib vs Vehicle within ", mes)
        dat_mes <- subset(dat, MES == mes)
        data_subset_for_tests <- data_for_tests[data_for_tests$MES == mes,]
        mes_stim_result <- perform_wilcox_test(data_subset_for_tests, as.formula('score ~ stim'))
        generate_violin_plot(dat_mes, feature, comp, paste0(feature, "_Stim_", mes, ".png"), group = "stim", p_value = mes_stim_result$p_value, custom_colors = c("#FF0000", "#0000FF"))
      results_df <- rbind(results_df, data.frame(Feature=feature, Comparison=comp, Group1="Selumetinib", Group2="Vehicle",
                                                  Median1=mes_stim_result$medians[1], Median2=mes_stim_result$medians[2], P_Value=mes_stim_result$p_value))
      }
  }
  # Save results to CSV
  write.csv(results_df, paste0(name, "_all_comparisons.csv"), row.names = FALSE)
}

# Call the function with Cdkn2a
feature_dir <- "neftel/MES/feature"
ifelse(!dir.exists(file.path(main_dir, feature_dir)), 
       dir.create(file.path(main_dir, feature_dir)), FALSE)
setwd(file.path(main_dir, feature_dir))
features <- c("Cdkn2a")
perform_wilcox_and_plots(features, tumor, "Cdkn2a")

# Call the function with G2.M.Score and G1.S.Score
feature_dir <- "neftel/MES/cellcyle"
ifelse(!dir.exists(file.path(main_dir, feature_dir)), 
       dir.create(file.path(main_dir, feature_dir)), FALSE)
setwd(file.path(main_dir, feature_dir))
features <- c("G2.M.Score", "G1.S.Score")
perform_wilcox_and_plots(features, tumor, "cellcycle")

#genesets
library(readxl)
library(nichenetr)
mek_activation <- read_excel("ref/Dry-Pratilas-Signature.xlsx")
mek_activation <- as.vector(mek_activation$gene) %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)]

mek_sens <- c("BRAF", "SHOC2", "SOS1", "SKP2") %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)]

mek_res <- c("ADAT2", "AMER1", "ANKMY2", "ARCN1", "ARL2", "ASNS", "ATP1A1", "ATR", "BCAS2",
           "BDP1", "BNIP1", "CCDC86", "CCT5", "CCT6A", "CDC27", "CDK1", "CEBPZ", "CENPE", 
           "CENPT", "CENPW", "CPSF4", "CSE1L", "CUL1", "DBR1", "DDB1", "DDX10", "DDX20", 
           "DDX21", "DDX23", "DDX24", "DDX41", "DDX52", "DHX37", "DHX8", "DKC1", "DNAJA3", 
           "DNTTIP2", "DSN1", "DYNC1I2", "DYNLRB1", "EFTUD1", "EIF2B5", "EIF3A", "EIF3D", 
           "EIF3E", "EIF3G", "EIF3I", "EIF6", "EPRS", "ERCC2", "EXOSC4", "EXOSC9", "FIP1L1", 
           "GCN1L1", "GLMN", "GMPS", "GNL2", "GPN3", "HEATR1", "INTS8", "ISG20L2", "KRI1", 
           "LAS1L", "MAK16", "MASTL", "MBTPS2", "MDN1", "MED1", "MED12", "MED7", "MPHOSPH6", 
           "MRPL35", "MRPS14", "MTOR", "MYBL2", "NABP2", "NAE1", "NDNL2", "NOL6", "NSMCE4A", 
           "NUBP1", "NUDCD3", "NUP54", "NUP98", "OPA1", "PAF1", "PAM16", "PCF11", "PELP1", 
           "POLR1E", "POLR2A", "POLR2H", "POLR3D", "POMP", "PPA1", "PRPF19", "PRPF31", "RABGGTA", 
           "RAD21", "RARS", "RASA2", "RBBP8", "RCL1", "RIOK1", "RIOK2", "RPAP3", "RPL17", 
           "RPL23", "RPL24", "RPL26", "RPL30", "RPL31", "RPL39", "RPL4", "RPL7", "RPS20", 
           "RPS24", "RPS27", "RPS6", "RRP12", "RRP9", "RRS1", "RUVBL2", "SAE1", "SAMM50", 
           "SKA1", "SMC2", "SMYD4", "SNAP23", "SNRNP200", "SOCS5", "SRP9", "SRPRB", "STAG2", 
           "SUSD1", "TACC3", "TAF1", "TAF2", "TAF6", "TAMM41", "TANGO6", "TARS2", "TCERG1", 
           "TEX10", "THOC3", "TMEM242", "TNPO3", "TTC4", "TUBGCP2", "TUBGCP3", "TUBGCP4", 
           "UFM1", "UQCRB", "UQCRFS1", "URI1", "UTP20", "UTP23", "UTP3", "VWA9", "WDR12", 
           "WDR18", "WDR5", "WDR55", "WDR82", "ZNHIT6") %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)]

# Define the gene lists and their names
gene_lists <- list(
  mek_activation = mek_activation,
  mek_sens = mek_sens,
  mek_res = mek_res)

  for (name in names(gene_lists)) {
   # Add Module Score
  tumor <- AddModuleScore(tumor, features = list(gene_lists[[name]]), name = name)
  }

# Generate the feature names from gene_lists keys
feature_names <- paste0(names(gene_lists), "1")

# Call the function with dynamically generated feature names
feature_dir <- "neftel/MES/genesets"
ifelse(!dir.exists(file.path(main_dir, feature_dir)), 
       dir.create(file.path(main_dir, feature_dir)), FALSE)
setwd(file.path(main_dir, feature_dir))
perform_wilcox_and_plots(feature_names, tumor, "genesets")

library(ggplot2)
library(dplyr)

# Create a summary data frame that calculates the proportion of Cycling cells by MES status for each sample
composition_data <- tumor@meta.data %>%
  filter(Cell_Cycle == "Cycling") %>% # Filter to include only cycling cells
  group_by(sample, MES) %>% # Group by sample and MES status
  summarise(Cycling_Count = n()) %>% # Count the number of cycling cells
  group_by(sample) %>% # Group by sample again
  mutate(Proportion_Cycling = Cycling_Count / sum(Cycling_Count)) # Calculate the proportion of cycling cells

# Plot the box plot showing the proportion of cycling cells in MES vs nonMES conditions
custom_colors <- c("#387ee0", "gray") # Define custom colors
boxplot <- ggplot(composition_data, aes(x = MES, y = Proportion_Cycling, fill = MES)) +
  geom_boxplot(width = 0.6, outlier.shape = 21, outlier.fill = "red", outlier.size = 2) +
  scale_fill_manual(values = custom_colors) +
  geom_jitter(width = 0.2, size = 1.5, color = "black", alpha = 0.5) + # Add jittered points
  labs(title = "Proportion of Cycling Cells in MES vs nonMES Conditions",
       x = "MES Condition",
       y = "Proportion of Cycling Cells") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16), # Center and bold title
    axis.text = element_text(size = 12), # Increase axis text size
    axis.title = element_text(size = 14, face = "bold"), # Bold axis titles
    axis.line = element_line(color = "black"), # Add axis lines
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_blank(), # Keep panel background minimal but retain axes
    plot.background = element_blank() # Remove plot background
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) # Show y-axis as percentage

ggsave("cycling_counts.pdf", boxplot, width = 7, height = 7)

#save rds
# saveRDS(tumor, file = paste0(main_dir, "/rds/tumor_neftel.rds"))
# tumor<-readRDS(paste0(main_dir, "/rds/tumor_neftel.rds"))

##################################
### ---------- dotplot  ---------- ###
##################################
feature_dir <- "neftel/MES/dotplot"
ifelse(!dir.exists(file.path(main_dir, feature_dir)), 
       dir.create(file.path(main_dir, feature_dir)), FALSE)
setwd(file.path(main_dir, feature_dir))
 
features <- sort(c("Fos", "Jun", "Junb", "Jund", "Nras","Dusp1"))

#by MES, split by MES
Idents(dat)= "MES"
dotplot <- DotPlot(dat, features = features, split.by = "stim", cols=c("#0000FF","#FF0000")) + 
              RotatedAxis()+ guides(color = guide_colorbar(title = 'Average Expression'))
ggsave(filename = "ERK_dotplot_MES_stim.png", plot = dotplot, width = 9, height = 5)

features <- sort(c("Ptprz1", "Olig1", "Olig2", "Pdgfra", "Sox10"))

#by MES, split by MES
Idents(dat)= "MES"
dotplot <- DotPlot(dat, features = features, split.by = "stim", cols=c("#0000FF","#FF0000")) + 
              RotatedAxis()+ guides(color = guide_colorbar(title = 'Average Expression'))
ggsave(filename = "De-differentiation_dotplot_MES_stim.png", plot = dotplot, width = 9, height = 5)
