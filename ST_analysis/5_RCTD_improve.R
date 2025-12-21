#De-convolute spot-level spatial data with scRNA-seq reference
library(spacexr)
library(Seurat)
library(ggplot2)
#Set up reference
ref <- readRDS("/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/nodomap_cleaned_20251025.RDS")
DefaultAssay(ref) <- "originalexp"
Idents(ref) <- "cluster_number"

#Extract counts/cluster information for RCTD reference
counts <- ref[["originalexp"]]@counts
##Change the rowname of reference counts (ENSEMBL ID) into gene symbol
mouse_gene_names<-read.delim("/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/mouse_genes_names_Ens100_rmdup.txt",row.names=3)
rownames(counts) <- rownames(mouse_gene_names[match(rownames(counts), mouse_gene_names$Gene.stable.ID),])

cluster <- as.factor(ref$cluster_number)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_originalexp
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

rm(ref)
gc()

#Set up query with spatial data
spatial <- readRDS("/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/nodose_st_RCTD_2401.RDS")
#remove old metadata 
spatial@meta.data$cell_types <- NULL
spatial <- SetIdent(spatial, value = "seurat_clusters")
spatial@meta.data <- spatial@meta.data[, !colnames(spatial@meta.data) %in% c("Sodium_channel", "Fibre_type", "Sensor_type", "Organ_projection")]
spatial@meta.data <- spatial@meta.data[, !colnames(spatial@meta.data) %in% c("first_type", "second_type", "spot_class", "first_class", "second_class", "min_score", "singlet_score", "conv_all", "conv_doublet")]
DefaultAssay(spatial) <- "RNA"
counts <- spatial[["RNA"]]@counts
##Get pixel data of integrated images 
# consider to offset coordinates so you can visualise all on the same plot 
coords <- list()
st.image <- Images(spatial)

for (i in seq_along(st.image)) {
  img_name <- st.image[i]
  tmp <- GetTissueCoordinates(spatial, image = img_name)
  
  # offset each slice along x-axis (or y if you prefer)
  offset_x <- (i - 1) * 7000  # shift each by 7000 pixels
  tmp$x <- tmp$x + offset_x
  
  coords[[img_name]] <- tmp
}
# combine all
coords <- do.call(rbind, coords)
colnames(coords) <- c("x", "y")
coords <- coords[, c("x", "y")]
rownames(coords) <- gsub(".*\\.", "", rownames(coords))

query <- SpatialRNA(coords, counts, colSums(counts))
print(dim(query@counts))
hist(log(query@nUMI,2)) # histogram of log2 nUMI per bead
print(head(query@coords))
barcodes <- colnames(query@counts) # pixels to be used (a list of barcode names). 
#plot_puck_continuous(query, barcodes, query@nUMI, ylimit = c(0,round(quantile(query@nUMI,0.9))), 
#                     title ='plot of nUMI')  

#Annotating and add cluster labels to the query (Spatial data)
RCTD <- create.RCTD(
  spatialRNA = query,
  reference = reference,
  max_cores = 8,
  gene_cutoff = 0.0001, # changes from 0.000125 !! min norm expression a gene must have to be kept for platform effect normalisation 
  fc_cutoff = 0.25, # changed from 0.5!!  min log fold change across cell types 
  gene_cutoff_reg = 1.5e-04, # chagned from 2e-04 !! min norm expression for main RCTD step
  fc_cutoff_reg = 0.5, # changed from 0.75 !! min log fold change for main RCTD step
  UMI_min_sigma = 100, # changedfrom 300!! when estimating the extra noise parameter only consider pixels with at least this many UMI so the estimate is not dominated by ultra noisy tiny pixels
)
#plot(RCTD@internal_vars$platform_effects)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
## update from here!! 
saveRDS(RCTD, "/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/RCTD_251025_paramrelax.RDS")

#doublet mode results
spatial <- AddMetaData(spatial, metadata = RCTD@results$results_df)
saveRDS(spatial, "/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/nodose_st_RCTD_paramrelax_251025.RDS")
Idents(spatial) <- "first_type"
#SpatialDimPlot(spatial, stroke = NA)
cat("Spot class overview:\n")
print(table(spatial$spot_class))

#full mode results
results <- RCTD@results
norm_weights <- normalize_weights(results$weights)
cell_type_name <- RCTD@cell_type_info$info[[2]]
spatialRNA <- RCTD@spatialRNA
setwd(dir = '/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/')
resultsdir <- 'RCTD_Plots_relaxed_params'
dir.create(resultsdir)
plot_weights(cell_type_name, spatialRNA, resultsdir, norm_weights)
plot_cond_occur(cell_type_name, resultsdir, norm_weights, spatialRNA)
plot_weights_doublet(cell_type_name, spatialRNA, resultsdir, results$weights_doublet,results$results_df)
