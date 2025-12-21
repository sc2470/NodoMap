## Trying a new method for neighbourhood analysis
library(igraph)
library(Seurat)
library(patchwork)
library(FNN)
library(ggplot2)
library(pheatmap)
library(MetBrewer)

## load data 
data_dir <- '/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/'
spatial<-readRDS(paste0(data_dir, 'nodose_st_RCTD_paramrelax_251025.RDS'))

# split into slices
spatial_list<-SplitObject(spatial, split.by = "orig.ident")
k <- 15
low_density_quantile <- 0.1 # filter bottom 10% by total neighbor weight. 

first_all  <- as.character(spatial$first_type)
second_all <- as.character(spatial$second_type)

global_annotations <- sort(unique(c(first_all, second_all)))
global_annotations <- global_annotations[global_annotations != "" & !is.na(global_annotations)]

global_annotations

spatial_list <- lapply(spatial_list, function(spatial) {
  
  # 1. remove reject cells 
  spatial <- subset(spatial, subset = spot_class != "reject")
  # 2. Get coordinates and prepare kNN
  coords <- GetTissueCoordinates(spatial)[, c("x", "y")]
  rownames(coords) <- rownames(GetTissueCoordinates(spatial))  # cell barcodes
  
  nn <- get.knn(coords, k = k)
  dists <- nn$nn.dist
  dists[dists == 0] <- 1e-6
  weights <- 1 / dists
  
  spatial@misc$knn_indices <- nn$nn.index
  spatial@misc$knn_weights <- weights
  
  # 3. compute inverse distance weighted counts 
  spot_class <- spatial$spot_class[rownames(coords)]
  first_type <- spatial$first_type[rownames(coords)]
  second_type <- spatial$second_type[rownames(coords)]
  
  stopifnot(length(spot_class) == nrow(coords))
  stopifnot(length(first_type) == nrow(coords))
  stopifnot(length(second_type) == nrow(coords))

  n_cells <- nrow(coords)
  annotations <- global_annotations
  weighted_counts <- matrix(0, nrow = n_cells, ncol = length(annotations))
  colnames(weighted_counts) <- annotations
  rownames(weighted_counts) <- rownames(coords)
  
  for (cell_i in 1:n_cells) {
    neighbors <- nn$nn.index[cell_i, ]
    neighbor_weights <- weights[cell_i, ]
    
    neighbor_classes <- spot_class[neighbors]
    neighbor_first <- first_type[neighbors]
    neighbor_second <- second_type[neighbors]
    
    for (j in 1:k) {
      if (neighbor_classes[j] == "singlet") {
        weighted_counts[cell_i, neighbor_first[j]] <- weighted_counts[cell_i, neighbor_first[j]] + neighbor_weights[j]
      } else if (neighbor_classes[j] %in% c("doublet_certain", "doublet_uncertain")) {
        weighted_counts[cell_i, neighbor_first[j]] <- weighted_counts[cell_i, neighbor_first[j]] + neighbor_weights[j]
        weighted_counts[cell_i, neighbor_second[j]] <- weighted_counts[cell_i, neighbor_second[j]] + neighbor_weights[j]
      }
    }
  }
  
  # 4. normalise weights per cell
  #weighted_counts_norm <- weighted_counts / rowSums(weighted_counts)  # composition-based
  #weighted_counts_norm[is.na(weighted_counts_norm)] <- 0  # handle rows with 0 total
  
  spatial[["neighbor_weights"]] <- CreateAssayObject(counts = t(weighted_counts))
  
  # 5. Plot total neighbor weight dist 
  total_nn <- rowSums(weighted_counts)
  df_plot <- data.frame(
    total_nn = total_nn
  )
  ggplot(df_plot, aes(x = total_nn)) +
    geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black", alpha = 0.8) +
    theme_minimal() +
    xlab("Total inverse-distance weighted neighbors") +
    ylab("Number of cells") +
    ggtitle(paste("Neighbor weight distribution -", spatial@project.name))
  
  return(spatial)
})

combined <- merge(spatial_list[[1]], y = spatial_list[-1], 
                  add.cell.ids = paste0("Slice", 1:length(spatial_list)), 
                  project = "SpatialCombined")

# scale, pca and cluster - without scaling 
DefaultAssay(combined) <- "neighbor_weights"

#combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
#combined@assays$neighbor_weights@scale.data <- as.matrix(combined[["neighbor_weights"]]@counts)
combined<-ScaleData(combined, features = rownames(combined), vars.to.regress= 'orig.ident')
combined <- RunPCA(combined, features = rownames(combined))
ElbowPlot(combined, ndims = 40)  # check variance explained

combined <- FindNeighbors(combined, dims = 1:10)  # adjust dims based on scree plot
combined <- FindClusters(combined, resolution = 0.2, cluster.name = 'neighbors' )  # adjust resolution to control cluster granularity

SpatialDimPlot(combined, group.by = 'neighbors')

neighbor_markers<-FindAllMarkers(combined, group.by = 'neighbors', only.pos = T)

######### Neighbourhood analysis plots ########
neighbourhood_names<-c("0" = "NGN/Glia",
                       "1" = "JGN",
                       "2" = "NGN",
                       "3" = "non-neuronal")
Idents(combined)<-'neighbors'
clusters <- as.character(Idents(combined))
neighborhood_label_vec <- neighbourhood_names[clusters]
names(neighborhood_label_vec) <- colnames(combined)
combined$neighborhood_label <- neighborhood_label_vec

## add back to main object 
## problems with differnt cell barcodes need to append slice to spatial barcodes 
labels_combined <- combined$neighborhood_label
labels_full <- rep("None", ncol(spatial))
names(labels_full) <- colnames(spatial)
labels_full[names(labels_combined)] <- labels_combined
spatial$neighborhood_label <- labels_full


custom_colors <- c(
  "NGN/Glia" = "#5D3A9B",  # deep purple
  "JGN" = "#E69F00",  # orange
  "NGN" = "#CC79A7",  # pink
  "non-neuronal" = "#0072B2",  # blue
  "None" = rgb(0.9,0.9,0.9,0.1)  # More transparent light grey (0.25 transparency)
)

plot<-SpatialDimPlot(combined, group.by = 'neighborhood_label', pt.size.factor = 1.5, cols = custom_colors, stroke = F, images = 'F1R1')
ggsave(filename = paste0(data_dir,"F1R1_neighbourhoods.png"),
       plot = plot, "png",dpi=1200,width=300,height = 300,units="mm")

plot<-SpatialDimPlot(combined, group.by = 'neighborhood_label', pt.size.factor = 1.5, cols = custom_colors, stroke = F, images = 'F2R2')
ggsave(filename = paste0(data_dir,"F2R2_neighbourhoods.png"),
       plot = plot, "png",dpi=1200,width=300,height = 300,units="mm")

plot<-SpatialDimPlot(combined, group.by = 'neighborhood_label', pt.size.factor = 1.5, cols = custom_colors, stroke = F, images = 'F3R3')
ggsave(filename = paste0(data_dir,"F3R3_neighbourhoods.png"),
       plot = plot, "png",dpi=1200,width=300,height = 300,units="mm")


#Heatmap to add 
avex<-round(log1p(data.frame(AverageExpression(combined, group.by = 'neighborhood_label')$neighbor_weights)), digits = 3)
zero_sum_rows <- rownames(avex)[rowSums(avex, na.rm = TRUE) == 0]
avex_filtered <- avex[rowSums(avex, na.rm = TRUE) != 0, ]
pheatmap(avex_filtered, 
         color = met.brewer('Hokusai2', 100, type = 'continuous'), 
         cluster_cols = T, 
         cluster_rows = T, 
         cellwidth = 8, 
         cellheight = 8, 
         treeheight_col = 0,
         scale='row',
         filname = paste0(data_dir, 'neighbourhood_heatmap.png'))



## testing plots 
coords_xy <- GetTissueCoordinates(spatial)[, c("x", "y")]

avg_dist <- rowMeans(nn$nn.dist)

library(ggplot2)
df <- data.frame(
  x = coords$x,
  y = coords$y,
  avg_dist = avg_dist
)

ggplot(df, aes(x = x, y = y, color = avg_dist)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  theme_minimal() +
  ggtitle("Average neighbor distance per cell")

plot(coords, pch = 19, col = "blue", asp = 1)

for (i in 1:nrow(nn$nn.index)) {
  for (j in 1:ncol(nn$nn.index)) {
    neighbor <- nn$nn.index[i, j]
    segments(coords[i, 1], coords[i, 2],
             coords[neighbor, 1], coords[neighbor, 2],
             col = rgb(0,0,0,0.2))  # semi-transparent edges
  }
}


# more plotting 

k <- ncol(nn$nn.index)
n_cells <- nrow(nn$nn.index)
annotations <- unique(cell_annot)

# Initialize matrix: cells × annotations
weighted_counts <- matrix(0, nrow = n_cells, ncol = length(annotations))
colnames(weighted_counts) <- annotations

# Compute inverse-distance weights
inv_dists <- 1 / nn$nn.dist
inv_dists[nn$nn.dist == 0] <- 1e6  # avoid division by zero

# Aggregate weights per annotation
for (i in 1:n_cells) {
  neighbors <- nn$nn.index[i, ]
  neighbor_annots <- cell_annot[neighbors]
  neighbor_weights <- inv_dists[i, ]
  
  weighted_counts[i, ] <- tapply(neighbor_weights, neighbor_annots, sum, default = 0)
}















