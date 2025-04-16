#Nodose spatial dimplot
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(scCustomize)
library(scales)

setwd("../../Downloads/2024_March_Nodomap/Spatial_Curio")
spatial <- readRDS("nodose_st_RCTD_2401.RDS")
DefaultAssay(spatial) <- "RNA"
levels(Idents(spatial))

#Set colour scale for different celltypes
celltype_col <- c("#FD6F87","#FC717F","#FA7476","#F8766D","#F47A5F","#F17E4F",
                  "#ED813C","#DE8C00","#D89000","#A3A500","#99A800","#8EAB00",
                  "#82AD00","#75AF00","#65B200","#53B400","#39B600","#00B70C",
                  "#00B931","#00BB45","#00BC56","#00BD64","#00BE71","#00BF7D",
                  "#00C089","#00C094","#00B3F1","#00B0F6","#00ACFB","#00A8FF",
                  "#06A4FF","#8893FF","#9590FF","#A18CFF","#AC88FF","#B584FF",
                  "#BF80FF","#C77CFF","#CF78FF","#D675FE","#DC71FA","#E26EF7",
                  "#E76BF3","#EC69EF","#F066EA","#F464E5","#F763E0","#FA62DB",
                  "#FC61D5","#FE61CF","#FF61C9","#FF61C2","#FF62BC")

celltypes <- c("EC1",
               "EC2",
               "EC3",
               "FB1",
               "FB2",
               "FB3",
               "FB4",
               "HC1",
               "HC2",
               "GC1",
               "GC2",
               "GC3",
               "MGC1",
               "MGC2",
               "MGC3",
               "MGC4",
               "MGC5",
               "MGC6",
               "SGC1",
               "SGC2",
               "SGC3",
               "SGC4",
               "SGC5",
               "SGC6",
               "SGC7",
               "SGC8",
               "JGN1",
               "JGN2",
               "JGN3",
               "JGN4",
               "JGN5",
               "NGN1",
               "NGN2",
               "NGN3",
               "NGN4",
               "NGN5",
               "NGN6",
               "NGN7",
               "NGN8",
               "NGN9",
               "NGN10",
               "NGN11",
               "NGN12",
               "NGN13",
               "NGN14",
               "NGN15",
               "NGN16",
               "NGN17",
               "NGN18",
               "NGN19",
               "NGN20",
               "NGN21",
               "NGN22")
Idents(spatial) <- factor(Idents(spatial), levels = celltypes)
levels(Idents(spatial))
names(celltype_col) <- levels(spatial)

x<-spatial@meta.data[,c(17,26)]
x$cell_type
x$cell_types[x$spot_class != 'singlet'] <- NA
colnames(x)<-c('spot_class', 'singlet_cells')
spatial <- AddMetaData(spatial, metadata = x$singlet_cells, col.name = "singlet_cells")

spatial@images$F1R1@coordinates[,1] <- spatial@images$F1R1@coordinates[,1]*-1
spatial@images$F1R1@coordinates[,2] <- spatial@images$F1R1@coordinates[,2]*-1

pdf("F1R1.pdf", width = 9.6, height = 9.6)
SpatialDimPlot(spatial, group.by = "singlet_cells", images = "F1R1", stroke = NA, 
               pt.size.factor = 1, image.alpha = 0) + NoLegend() + 
  scale_fill_manual(values = celltype_col, na.value="#d6d6d6")
dev.off()

pdf("F2R2.pdf", width = 9.6, height = 9.6)
SpatialDimPlot(spatial, group.by = "singlet_cells", images = "F2R2", stroke = NA, 
               pt.size.factor = 1, image.alpha = 0) + NoLegend() + 
  scale_fill_manual(values = celltype_col, na.value="#d6d6d6")
dev.off()

pdf("F3R3.pdf", width = 9.6, height = 9.6)
SpatialDimPlot(spatial, group.by = "singlet_cells", images = "F3R3", stroke = NA, 
               pt.size.factor = 1, image.alpha = 0) + NoLegend() + 
  scale_fill_manual(values = celltype_col, na.value="#d6d6d6")
dev.off()








