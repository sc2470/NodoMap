# ============================================
# Script: 2_subset_and_save.R
# Purpose: Remove NGN5 update metadata
# ============================================

# Load required packages ----
library(Seurat)
library(ggplot2)
library(stringr)

# Set up paths and config ----
data_dir <- '/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/'

# Load data ----
data <- readRDS(paste0(data_dir, 'nodose-nucseq.RDS'))
diet <- readRDS(paste0(data_dir, 'diet_nodomap_metadata_updated.RDS'))

# Analysis ----
# # Step 0: move diet metadata to full object 
diet_meta <- diet@meta.data
data <- AddMetaData(data, diet_meta)
DimPlot(data, group.by = 'cell_types')
rm(diet)

# Step 1: remove LQ cluster
subset_data <- subset(data, subset = cell_types != 'NGN5')
 
# Step 2: Redraw umap
DimPlot(subset_data, group.by = 'cell_types')

# Step 3: Rename nodose clusters
cn <- as.character(subset_data$Seurat_cluster)
nums <- suppressWarnings(as.integer(str_match(cn, "^NGN(\\d+)")[,2]))
nums_new <- ifelse(!is.na(nums) & nums >= 6, nums - 1, nums)
prefix_new <- ifelse(!is.na(nums_new), paste0("NGN", nums_new), NA_character_)
suffix <- ifelse(!is.na(nums),
                 str_replace(cn, "^NGN\\d+", ""), 
                 NA_character_)
cn_new <- ifelse(!is.na(nums),
                 paste0(prefix_new, suffix),
                 cn)

subset_data$cluster_number <- factor(cn_new, levels = sort(unique(cn_new)))
Idents(subset_data) <- subset_data$cluster_number

# Step 4: clean up metadata
subset_data@meta.data <- subset_data@meta.data[, !grepl("^integrated_snn", colnames(subset_data@meta.data))]
subset_data@meta.data <- subset_data@meta.data[, !colnames(subset_data@meta.data) %in% c("cell_types", "CELLEX_markers", "Seurat_cluster")]
subset_data@meta.data <- subset_data@meta.data[, !colnames(subset_data@meta.data) %in% c("Organ_projection", "Sodium_channel", "Fibre_type", "Sensor_type", "NPs", "NTs", "cluster_names")]
subset_data$seurat_clusters <- droplevels(subset_data$seurat_clusters)

# Save output ----
saveRDS(subset_data, file.path(data_dir, "nodomap_cleaned_20251025.RDS"))
