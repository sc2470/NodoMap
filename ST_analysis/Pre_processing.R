#Filter the aligned Curio nodose spatial transcriptomic data
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#load the data
##F1R1
F1R1 <- readRDS("../nodose_st_230918/F1R1_seurat.rds")
DefaultAssay(F1R1) = "RNA"

summary(F1R1@meta.data$nCount_RNA)
SpatialFeaturePlot(F1R1, "nCount_RNA", stroke = 0)

#Filter the cells with low counts
F1R1 <- subset(F1R1, subset = nCount_RNA > 100)
SpatialFeaturePlot(F1R1, "nCount_RNA", stroke = 0)
summary(F1R1$nCount_RNA)

#Filter the cells with low features
F1R1 <- subset(F1R1, subset = nFeature_RNA > 190)
SpatialFeaturePlot(F1R1, "nFeature_RNA", stroke = 0)
summary(F1R1$nFeature_RNA)

F1R1 <-AddMetaData(F1R1, "LNG", col.name = "Position")
names(F1R1@images) <- "F1R1"
F1R1 <- NormalizeData(F1R1)
saveRDS(F1R1, "F1R1_processed_20230921.rds")

##F2R2
F2R2 <- readRDS("../nodose_st_230918/F2R2_seurat.rds")
DefaultAssay(F2R2) = "RNA"

summary(F2R2@meta.data$nCount_RNA)
SpatialFeaturePlot(F2R2, "nCount_RNA", stroke = 0)

F2R2 <- subset(F2R2, subset = nCount_RNA >100)
SpatialFeaturePlot(F2R2, "nCount_RNA", stroke = 0)
summary(F2R2$nCount_RNA)

F2R2 <- subset(F2R2, subset = nFeature_RNA > 190)
SpatialFeaturePlot(F2R2, "nFeature_RNA", stroke = 0)
summary(F2R2$nFeature_RNA)


F2R2 <-AddMetaData(F2R2, "LNG", col.name = "Position")
names(F2R2@images) <- "F2R2"
F2R2 <- NormalizeData(F2R2)
saveRDS(F2R2, "F2R2_processed_20230921.rds")

##F3R3
F3R3 <- readRDS("../nodose_st_230918/F3R3_seurat.rds")
DefaultAssay(F3R3) = "RNA"

summary(F3R3$nCount_RNA)
SpatialFeaturePlot(F3R3, "nCount_RNA", stroke = 0)

F3R3 <- subset(F3R3, subset = nCount_RNA >100)
SpatialFeaturePlot(F3R3, "nCount_RNA", stroke = 0)
summary(F3R3$nCount_RNA)

#Filter the cells with low features
F3R3 <- subset(F3R3, subset = nFeature_RNA > 190)
SpatialFeaturePlot(F3R3, "nFeature_RNA", stroke = 0)
summary(F3R3$nFeature_RNA)

F3R3 <- AddMetaData(F3R3, "RNG", col.name = "Position")
names(F3R3@images) <- "F3R3"
F3R3 <- NormalizeData(F3R3)
saveRDS(F3R3, "F3R3_processed_20230921.rds")

