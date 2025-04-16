library(scater)
library(DropletUtils)
library(scDblFinder)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(hdf5r)
library(tidyr)

data<-list()
setwd("~/mouse-nuc-seq/nodose/data/processed")
files<-dir(pattern = ".RDS")
for (i in files)
{
  data[[i]]<-readRDS(i)
}

data<-lapply(X=data, assay = 'originalexp', FUN = SCTransform)
#saving the SCT datasets
saveRDS(data, "230401/SCT-datasets.RDS")

features <- SelectIntegrationFeatures(object.list = data, nfeatures = 3000)
data<-PrepSCTIntegration(object.list = data, anchor.features = features)
data.anchors<-FindIntegrationAnchors(data, normalization.method = "SCT", anchor.features = features)
saveRDS(data.anchors, "230401/anchors.RDS")
data.combined.sct<-IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
data.combined.sct<-RunPCA(data.combined.sct, verbose = FALSE)

ElbowPlot(data.combined.sct,ndims = 50)

data.combined.sct<-RunUMAP(data.combined.sct, reduction = "pca", dims = 1:30)

data.combined.sct<-FindNeighbors(data.combined.sct,reduction = "pca",dims = 1:30)
data.combined.sct<-FindClusters(data.combined.sct,resolution = 0.8)

DimPlot(data.combined.sct, reduction = "umap")
DimPlot(data.combined.sct, reduction = "umap", group.by = "Dataset")
DefaultAssay(data.combined.sct)

DefaultAssay(data.combined.sct)<-'originalexp'

FeaturePlot(data.combined.sct, "nCount_originalexp")


gene_names<-read.delim("../../../../mouse-nuc-seq/mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes<-function(gene_symbol){
  return(as.character(gene_names[gene_symbol,1]))
}



saveRDS(data.combined.sct, '230401/nodose-integrated.RDS')
