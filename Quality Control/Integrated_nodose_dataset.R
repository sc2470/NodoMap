library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

#load integrated nodose data
data <- readRDS("/home/rstudio/shared_data/nodose/latest-integration-230403/nodose-integrated.RDS")

DefaultAssay(data) = "originalexp"

#PC30, resolution = 0.8

pdf("figure_QC_Vln.pdf", width = 12, height = 6)

VlnPlot(data, "nCount_originalexp", log = T, pt.size = 0) + NoLegend()
VlnPlot(data, "nFeature_originalexp", log = T, pt.size = 0) + NoLegend()
VlnPlot(data, "subsets_Mito_percent", log = T, pt.size = 0) + NoLegend()

dev.off()

pdf("figure_QC_Feature.pdf", width = 6, height = 5)

FeaturePlot(data, "nCount_originalexp", label = T)
FeaturePlot(data, "nFeature_originalexp", label = T)
FeaturePlot(data, "subsets_Mito_percent", label = T)

dev.off()

pdf("figure_QC_Scatter.pdf", width = 6, height = 5)

FeatureScatter(data, "nCount_originalexp", "subsets_Mito_percent")
FeatureScatter(data, "nCount_originalexp", "nFeature_originalexp")

dev.off()

summary(data$nCount_originalexp)
summary(data$nFeature_originalexp)
summary(data$subsets_Mito_percent)

#umaps
pdf("figrue_umap.pdf", width = 10, height = 10)

DimPlot(data, reduction = "umap", label = T, raster = F, pt.size = 0.1) + NoLegend()

DimPlot(data, reduction = "umap", group.by = "Dataset", raster = F, pt.size = 0.1) + 
  scale_color_manual(
    labels = c("Bai et al., 2019, scSeq", 
               "Buchanan et al., 2022, scSeq", 
               "Cheng et al., 2023, nucSeq", 
               "Kupari et al., 2019, scSeq", 
               "Zhao et al., 2022, scSeq"), 
    values = c('#CC4E5C', '#DAA520', 'darkgreen', '#318CE7', '#fbefef'))

DimPlot(data, reduction = "umap", group.by = "Position", repel = T, raster = F, pt.size = 0.1,
        cols = c('left'='#31AABF', 'right'='#E566B2'), na.value = "grey99")


DimPlot(data, reduction = "umap", group.by = "Nutr.cond", repel = T, raster = F, pt.size = 0.1,
        cols = c('adlib'='#488B27', 'fast'='#8C00E5'), na.value = "grey99")

dev.off()

#Function to search by gene names
mouse_gene_names <- read.delim("/home/rstudio/mouse_genes_names_Ens100_rmdup.txt", row.names=3)
genes_to_ens <- function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}


#Cell types
library(RColorBrewer)
palette(brewer.pal(10,"Spectral"))

nonneurongenes <- c("Ptprc", "Lum", "Ebf2", "Pdgfra", "Cspg4", "Emcn", "Cldn5")
glialgenes <- c("Bcan", "Kcnj10", "Mbp", "Sox10", "Foxj1", "Aqp4")
neurongenes <- c("Fos", "Rbfox3", "Slc17a6", "Slc32a1", "Slc18a3")
NGN <- c("Phox2b", "Ntrk2", "P2rx2", "Eya1")
JGN <- c("Prdm12", "Prrxl1", "Gm13425")

pdf("figure_geneofinterests.pdf", width = 6, height = 5)

for(i in 1:length(nonneurongenes))
{ print(i)
  p <- FeaturePlot(data, genes_to_ens(nonneurongenes[i]), label = T, 
                   order = T,cols=c(rgb(0.9,0.9,0.9,0.1),palette()[i])) + ggtitle(nonneurongenes[i])
plot(p)
}

for(i in 1:length(glialgenes))
{ print(i)
  p<-FeaturePlot(data, genes_to_ens(glialgenes[i]), label = T, 
                 order = T,cols=c(rgb(0.9,0.9,0.9,0.1),palette()[i])) + ggtitle(glialgenes[i])
  plot(p)
}

for(i in 1:length(neurongenes))
{ print(i)
  p<-FeaturePlot(data, genes_to_ens(neurongenes[i]), label = T, 
                 order = T,cols=c(rgb(0.9,0.9,0.9,0.1),palette()[i])) + ggtitle(neurongenes[i])
  plot(p)
}

for(i in 1:length(NGN))
{ print(i)
  p<-FeaturePlot(data, genes_to_ens(NGN[i]), label = T,  
                 order = T,cols=c(rgb(0.9,0.9,0.9,0.1),palette()[i])) + ggtitle(NGN[i])
  plot(p)
}

for(i in 1:length(JGN))
{ print(i)
  p<-FeaturePlot(data, genes_to_ens(JGN[i]), label = T, 
                 order = T,cols=c(rgb(0.9,0.9,0.9,0.1),palette()[i])) + ggtitle(JGN[i])
  plot(p)
}

dev.off()

#Number of cells in each clusters/study/batch
tab <- table(data@active.ident)
write.table(tab, file = "total_cell_cluster.txt", quote = F, sep = "\t", col.names = F)

tab <- table(data@meta.data$Dataset)
write.table(tab, file = "total_cell_dataset.txt", quote = F, sep = "\t", col.names = F)

tab <- table(data@meta.data$Experiment)
write.table(tab, file = "total_cell_batch.txt", quote = F, sep = "\t", col.names = F)

#Number of cells in each cluster for each study
tab <- table(paste0(data@active.ident, data@meta.data$Dataset))
write.table(tab, file = "cell_cluster_dataset.txt", quote = F, sep = "\t", col.names = F)

#Number of cells in each cluster for each batch
tab <- table(paste0(data@active.ident, data@meta.data$Experiment))
write.table(tab, file = "cell_cluster_batch.txt", quote = F, sep = "\t", col.names = F)

#Number of left & right cells
tab <- table(paste(data@meta.data$Experiment, data@meta.data$Position))
write.table(tab, file = "cell_left_right.txt", quote = F, sep = "\t", col.names = F)

#Number of fed & fasted cells
tab <- table(paste(data@meta.data$Dataset, data@meta.data$Nutr.cond))  
write.table(tab, file = "cell_fed_fasted.txt", quote = F, sep = "\t", col.names = F)
