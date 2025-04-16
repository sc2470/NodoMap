library(Seurat)
library(spacexr)
library(ggplot2)

RCTD <- readRDS("RCTD_231031.RDS")
spatial <- readRDS("nodose_st_RCTD.RDS")
DefaultAssay(spatial) <- "RNA"
table(spatial$spot_class)

#Add cluster name to spatial data
spatial <- SetIdent(spatial, value = "first_type")
clustermarkers <- c("EC1",
                    "SGC1",
                    "SGC2",
                    "MGC1",
                    "SGC3",
                    "NGN1",
                    "HC1",
                    "SGC4",
                    "SGC5",
                    "FB1",
                    "SGC6",
                    "MGC2",
                    "FB2",
                    "NGN2",
                    "GC1",
                    "NGN3",
                    "MGC3",
                    "NGN4",
                    "SGC7",
                    "MGC4",
                    "NGN5",
                    "NGN6",
                    "SGC8",
                    "JGN1",
                    "NGN7",
                    "GC2",
                    "JGN2",
                    "FB3",
                    "NGN8",
                    "NGN9",
                    "EC2",
                    "NGN10",
                    "NGN11",
                    "NGN12",
                    "NGN13",
                    "EC3",
                    "NGN14",
                    "NGN15",
                    "JGN3",
                    "HC2",
                    "NGN16",
                    "JGN4",
                    "NGN17",
                    "NGN18",
                    "FB4",
                    "GC3",
                    "NGN19",
                    "JGN5",
                    "MGC5",
                    "MGC6",
                    "NGN20",
                    "NGN21",
                    "NGN22")

names(clustermarkers)<-levels(spatial)
spatial<-RenameIdents(spatial, clustermarkers)
spatial <- AddMetaData(object = spatial, metadata = Idents(spatial), col.name = "cell_types")
saveRDS(spatial,"nodose_st_RCTD_2401.RDS")

pdf("Spatial_dimplot_RCTD.pdf", width = 10, height = 10)
SpatialDimPlot(spatial,images = "F1R1",stroke = NA)
SpatialDimPlot(spatial,images = "F2R2",stroke = NA)
SpatialDimPlot(spatial,images = "F3R3",stroke = NA)
dev.off()

pdf("Dimplot_RCTD.pdf", width = 50, height = 25)
DimPlot(spatial,label = T,label.size = 6,repel = T,raster = F,pt.size = 3) + 
  guides(color = guide_legend(override.aes = list(size=3), ncol = 1)) +
  theme(axis.title = element_text(face = "bold", size = 30),
        axis.text = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(face = "bold", size = 30, margin = margin(t = 2, b = 2, unit = "pt"))) +
  coord_fixed()
dev.off()

spatial_true <- subset(spatial,subset = spot_class != "reject")
pdf("Dimplot_RCTD_true.pdf", width = 50, height = 25)
DimPlot(spatial_true,label = T,label.size = 6,repel = T,raster = F,pt.size = 3) + 
  guides(color = guide_legend(override.aes = list(size=3), ncol = 1)) +
  theme(axis.title = element_text(face = "bold", size = 30),
        axis.text = element_text(face = "bold", size = 30),
        axis.line = element_line(linewidth = 1),
        axis.ticks = element_line(linewidth = 1),
        legend.text = element_text(face = "bold", size = 30, margin = margin(t = 2, b = 2, unit = "pt"))) +
  coord_fixed()
dev.off()

pdf("Spatial_dimplot_RCTD_true.pdf", width = 10, height = 10)
SpatialDimPlot(spatial_true,images = "F1R1",stroke = NA)
SpatialDimPlot(spatial_true,images = "F2R2",stroke = NA)
SpatialDimPlot(spatial_true,images = "F3R3",stroke = NA)
dev.off()


DimPlot(spatial, reduction = "umap", group.by = "Position")






