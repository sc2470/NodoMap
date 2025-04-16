#Integrate nodose spatial transcriptomic dataset
library(Seurat)
library(patchwork)
library(dplyr)

#Load the aligned data
nodose.st <- list()
files <- list.files(path = ".", pattern = ".rds")
for (i in files) {
  nodose.st[[i]] <- readRDS(i)
}

#SC transform each .rds file
nodose.st <- lapply(nodose.st, assay = "RNA", FUN = SCTransform)

features <- SelectIntegrationFeatures(object.list = nodose.st, nfeatures = 3000)
nodose.st <- PrepSCTIntegration(object.list = nodose.st, anchor.features = features)

nodose.anchors <- FindIntegrationAnchors(object.list = nodose.st, normalization.method = "SCT", anchor.features = features)

saveRDS(nodose.anchors, "anchors_230921.RDS")

nodose.combined.sct <- IntegrateData(anchorset = nodose.anchors, normalization.method = "SCT")
nodose.combined.sct <- RunPCA(nodose.combined.sct, verbose=FALSE)

ElbowPlot(nodose.combined.sct, ndims = 50)

nodose.combined.sct <- RunUMAP(nodose.combined.sct, reduction = "pca", dims = 1:40)
nodose.combined.sct <- FindNeighbors(nodose.combined.sct, reduction = "pca", dims = 1:40)
nodose.combined.sct <- FindClusters(nodose.combined.sct, resolution = 0.2)

DimPlot(nodose.combined.sct, reduction = "umap", label = T)
DimPlot(nodose.combined.sct, reduction = "umap", group.by = "orig.ident")
DimPlot(nodose.combined.sct, reduction = "umap", group.by = "Position")
SpatialDimPlot(nodose.combined.sct, stroke = NA)

saveRDS(nodose.combined.sct, "nodose_st_combined_sct.RDS")













