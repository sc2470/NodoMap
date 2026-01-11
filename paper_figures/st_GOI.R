library(Seurat)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(stringr)

#load spatial seurat object
data <- readRDS("nodose_st_RCTD_251025_metadata_updated.RDS")
DefaultAssay(data) <- "RNA"
View(data@meta.data)

##Reorder levels by cell types
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
               "NGN21")
data <- SetIdent(data, value = "first_type")
levels(Idents(data))
Idents(data) <- factor(Idents(data), levels = celltypes)

colnames(data@meta.data)[26]
data$non_neuronal <- ifelse(
  data$singlet_cells %in% c("MGC2","MGC3","SGC1","EC1","NGN11","NGN3","EC2","FB4","NGN5","EC3","FB2"),
  as.character(data$singlet_cells),
  NA)
View(data@meta.data)
unique(data$non_neuronal)

#spatial dimplot with non_neuronal neighbourhood
celltype_col <- c("#FD6F87","#FC717F","#FA7476","#F8766D","#F47A5F","#F17E4F",
                  "#ED813C","#DE8C00","#D89000","#A3A500","#99A800","#8EAB00",
                  "#82AD00","#75AF00","#65B200","#53B400","#39B600","#00B70C",
                  "#00B931","#00BB45","#00BC56","#00BD64","#00BE71","#00BF7D",
                  "#00C089","#00C094","#00B3F1","#00B0F6","#00ACFB","#00A8FF",
                  "#06A4FF","#8893FF","#9590FF","#A18CFF","#AC88FF","#B584FF",
                  "#BF80FF","#C77CFF","#CF78FF","#D675FE","#DC71FA","#E26EF7",
                  "#E76BF3","#EC69EF","#F066EA","#F464E5","#F763E0","#FA62DB",
                  "#FC61D5","#FE61CF","#FF61C9","#FF61C2")

names(celltype_col)<-levels(data)

data$non_neuronal <- factor(data$non_neuronal, levels = celltypes)

pdf("RCTD_nn_dimplot.pdf", width = 5, height = 5)
SpatialDimPlot(data,group.by = "non_neuronal",images = "F2R2",
               stroke = NA,image.alpha = 0,pt.size.factor = 1.3) +
  scale_fill_manual(values = celltype_col, na.value = rgb(.95,.95,.95))
dev.off()

#Identify non-neuronal neighbourhood cluster markers
#Find markers of every clusters, only show positive ones
data.marker <- FindAllMarkers(data, only.pos = T, min.pct = 0.25)

pdf("RCTD_nn_featureplot.pdf", width = 9.6, height = 9.6)
SpatialFeaturePlot(data,"Fabp4",images = "F2R2",pt.size.factor = 2,
                   alpha = c(0.5,1),image.alpha = 0) + 
  scale_fill_gradient(high = "springgreen4",low = "#d6d6d6")

SpatialFeaturePlot(data,"Hbb-bt",images = "F2R2",pt.size.factor = 2,
                   alpha = c(0.5,1),image.alpha = 0) + 
  scale_fill_gradient(high = "darkorange3",low = "#d6d6d6")
dev.off()



# GOI spatial feature plot ------------------------------------------------
goi <- c("Cckar", "Cnr1", "Gabra1", "Htr3a", 
         "Cartpt", "Lypd1", "Nnat", "Nrsn1", "Olfm1", "Scg2", 
         "Slc25a3", "Slc25a39")
df <- as.data.frame(goi)
goi_col <- c("violetred3","mediumorchid4","orangered4","deeppink4",
             "chartreuse4","mediumseagreen","dodgerblue3","firebrick4","coral3","slateblue3",
             "steelblue3","olivedrab")

pdf("figure_st_goi.pdf", width = 8, height = 4)
for (i in 1:length(df$goi)) {
  p <- SpatialFeaturePlot(data, df$goi[i], keep.scale = "all",pt.size.factor = 2,
                          alpha = c(0.5,1),image.alpha = 0) & 
    scale_fill_gradient(high = goi_col[i],low = "#d6d6d6")
  plot(p)
}
dev.off()



















