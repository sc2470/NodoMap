#De-convolute spot-level spatial data with scRNA-seq reference
library(spacexr)
library(Seurat)
library(ggplot2)
#Set up reference
ref <- readRDS("../2023_April_Integration/nodose_integrated_PC30-60_res0.6-2/nodose-integrated-pcs-50-resolution-1.4.RDS")
ref <- UpdateSeuratObject(ref)
DefaultAssay(ref) <- "originalexp"
Idents(ref) <- "integrated_snn_res.1.4"

#Extract counts/cluster information for RCTD reference
counts <- ref[["originalexp"]]@counts
##Change the rowname of reference counts (ENSEMBL ID) into gene symbol
mouse_gene_names<-read.delim("../mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}
rownames(counts) <- rownames(mouse_gene_names[match(rownames(counts), mouse_gene_names$Gene.stable.ID),])

cluster <- as.factor(ref$seurat_clusters)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_originalexp
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

#Set up query with spatial data
spatial <- readRDS("nodose_st_combined_sct.RDS")
DefaultAssay(spatial) <- "RNA"
counts <- spatial[["RNA"]]@counts
##Get pixel data of integrated images
coords <- list()
st.image <- Images(spatial)
for (i in st.image) {
  coords[[i]] <- GetTissueCoordinates(spatial, image = i)
}
coords <- do.call(rbind,coords)
colnames(coords) <- c("x","y")
coords[is.na(colnames(coords))] <- NULL
rownames(coords) <- gsub(".*\\.","",rownames(coords))
query <- SpatialRNA(coords, counts, colSums(counts))

#Annotating and add cluster labels to the query (Spatial data)
RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
saveRDS(RCTD, "RCTD_231031.RDS")

#doublet mode resutls
spatial <- AddMetaData(spatial, metadata = RCTD@results$results_df)
saveRDS(spatial, "nodose_st_RCTD.RDS")
Idents(spatial) <- "first_type"
SpatialDimPlot(spatial, stroke = NA)

#full mode results
results <- RCTD@results
norm_weights <- normalize_weights(results$weights)
cell_type_name <- RCTD@cell_type_info$info[[2]]
spatialRNA <- RCTD@spatialRNA
resultsdir <- 'RCTD_Plots'
dir.create(resultsdir)
plot_weights(cell_type_name, spatialRNA, resultsdir, norm_weights)
plot_cond_occur(cell_type_name, resultsdir, norm_weights, spatialRNA)
plot_weights_doublet(cell_type_name, spatialRNA, resultsdir, results$weights_doublet,results$results_df)







