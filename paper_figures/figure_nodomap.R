library(Seurat)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(RColorBrewer)
library(viridis)

data<-readRDS("/home/rstudio/2023_April_Integration/nodose_integrated_PC30-60_res0.6-2/nodose-integrated-pcs-50-resolution-1.4.RDS")
DefaultAssay(data) = "originalexp"
data <- SetIdent(data, value = "integrated_snn_res.1.4")

mouse_gene_names<-read.delim("../../mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

#Figure 1: Dimplot
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

names(clustermarkers)<-levels(data)
data<-RenameIdents(data, clustermarkers)
data <- AddMetaData(object = data, metadata = Idents(data), col.name = "cell_types")

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
               "NGN21",
               "NGN22")
Idents(data) <- factor(Idents(data), levels = celltypes)

legend <- c("EC1.Cldn5/Ly6c1/Ly6a",
            "EC2.Selp/Fabp4/Vwf",
            "EC3.Podxl/Depp1/Ptprb",
            "FB1.Lum/Ccl11/Smoc2",
            "FB2.Acta2/Myh11/Rgs5",
            "FB3.Igfbp6/Thbs4/Islr",
            "FB4.Hbb-bs/Bnc2/Ebf2",
            "HC1.C1qc/C1qa/C1qb",
            "HC2.Cd52/Lsp1/Rac2",
            "GC1.Kcna1/Fxyd3/Mbp",
            "GC2.Sostdc1/Plekhb1/Tsc22d4",
            "GC3.C1qc/Tyrobp/Lyz2",
            "MGC1.Ncmap/Fam178b/Cldn19",
            "MGC2.Prx/Pllp/Bcas1",
            "MGC3.Slc36a2/Nr4a2/Ugt8a",
            "MGC4.Mt2/Tcim/Drp2",
            "MGC5.Smoc2/Ccl11/Lum",
            "MGC6.Ifitm1/Kcnj8/Rgs5",
            "SGC1.Mmd2/Ednrb/Fabp7",
            "SGC2.Cxcl10/Acsbg1/Phgdh",
            "SGC3.Bcan/Fbln2/Fbln5",
            "SGC4.Atf3/Emp1/Cebpd",
            "SGC5.Pou3f1/Col1a2/Gas7",
            "SGC6.Ifit3/Cdh19/Isg15",
            "SGC7.G0s2/Sdc4/Ccn1",
            "SGC8.Ttyh1/Ptprz1/Sfrp5",
            "JGN1.Mrgprd/Tmem45b/Grik1",
            "JGN2.Gfra3/Tac1/Tmem255a",
            "JGN3.Trappc3l/Nptx1/Tafa1",
            "JGN4.Wfdc2/C1ql4/Rarres1",
            "JGN5.Trpm8/Foxp2/Gsg1l",
            "NGN1.C1ql3/Adcyap1/St18",
            "NGN2.Htr3b/Htr3a/Hs3st4",
            "NGN3.Cysltr2/Chrnb3/Kcnip4",
            "NGN4.Uts2b/Cckar/Vip",
            "NGN5.Rtn1/Ngfr/Clu",
            "NGN6.Kcng1/Vmn1r85/Trpv1",
            "NGN7.Lypd1/Ddc/Tafa2",
            "NGN8.Rbp4/Gda/Kcnk9",
            "NGN9.Trpa1/Nos1/Efcab6",
            "NGN10.Miat/Snhg11/Meg3",
            "NGN11.Olfm3/Tmem233/P2ry1",
            "NGN12.Sprr1a/Ecel1/Nts",
            "NGN13.Slc18a3/Pappa2/Runx3",
            "NGN14.Lox/Pkib/Kcnk3",
            "NGN15.Slc17a7/Hapln4/Lgi3",
            "NGN16.Bmp3/Gal/Sost",
            "NGN17.Gabra1/Chodl/Gabrb2",
            "NGN18.Thsd7b/Cacng5/Gata3",
            "NGN19.Gabra1/Gabrb2/Chrm2",
            "NGN20.Slc6a2/Npy/Hand2",
            "NGN21.Amigo2/Glp1r/Npy2r",
            "NGN22.Olfr78/Runx3/Avpr1a")

names(legend)<-levels(data)
data<-RenameIdents(data,legend)
data <- AddMetaData(data, metadata = Idents(data), col.name = "Wilcox_markers")
saveRDS(data, "nodose_sc_integrated_Wicox_2401.rds")

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
PalettePlot(palette = celltype_col)

#Dimplot
x1 <- DimPlot_scCustom(data, colors_use = celltype_col, reduction = "umap", raster = FALSE, label = F) 
x2 <- x1+NoLegend()+
  theme(axis.title = element_blank(),
        axis.text = element_blank()) + coord_fixed()

pdf("figure_dimplot.pdf", width = 10, height = 10)
print(x2)
dev.off()

x3 <- x1+NoAxes()+
  guides(color = guide_legend(override.aes = list(size=7), ncol = 5)) +
  theme(legend.text = element_text(size = 24, margin = margin(t = 3, b = 3, unit = "pt"), color = "grey40", family = 'Helvetica'))

pdf("figure_dimplot_legend.pdf", width = 26.5, height = 5)
print(x3)
dev.off()

#NGN
summary(GetAssayData(data, slot = "data")[genes_to_ens("Phox2b"),])
pdf("figure_Phox2b.pdf", width = 3, height = 3)
FeaturePlot(data, genes_to_ens("Phox2b"), raster = F) + NoAxes() + ggtitle(NULL) + 
  scale_colour_gradient(high = "violetred4", low = rgb(0.9,0.9,0.9,0.1), limits=c(0,6)) +
  theme(legend.text = element_text(size = 10, family = "Helvetica"),legend.spacing = unit(1.0, 'cm')) + 
  coord_fixed()
dev.off()

#Heatmap
NGN <- c("NGN1",
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

Idents(data) <- "cell_types"
neurons <- subset(data, idents = NGN)
DefaultAssay(neurons) = "originalexp"
levels(Idents(neurons))

all.genes <- rownames(neurons)
neurons <- ScaleData(neurons, features = all.genes)

#Average expression for each identity class
NGN <- AverageExpression(neurons,assays = "originalexp",group.by = "ident",return.seurat = T)
#Pull the list of colors used in DimPlot
nodose_col <- c("#8893FF","#9590FF","#A18CFF","#AC88FF","#B584FF","#BF80FF",
                "#C77CFF","#CF78FF","#D675FE","#DC71FA","#E26EF7","#E76BF3",
                "#EC69EF","#F066EA","#F464E5","#F763E0","#FA62DB","#FC61D5",
                "#FE61CF","#FF61C9","#FF61C2","#FF62BC")
PalettePlot(palette = nodose_col)

markers <- c("Adcyap1",
             "Htr3a",
             "Cysltr2",
             "Cckar",
             "Clu",
             "Kcng1",
             "Lypd1",
             "Rbp4",
             "Trpa1",
             "Miat",
             "Tmem233",
             "Sprr1a",
             "Slc18a3",
             "Pkib",
             "Slc17a7",
             "Sost",
             "Gabra1",
             "Thsd7b",
             "Chrm2",
             "Slc6a2",
             "Amigo2",
             "Olfr78")

pdf("figure_heatmap.pdf", width = 11, height = 12.5)
DoHeatmap(NGN, features = genes_to_ens(markers), size = 10, raster = F, 
          angle = 90, draw.lines = 0, group.colors = nodose_col) +
  scale_y_discrete(labels = rev(markers)) + 
  scale_fill_viridis()+
  guides(color= "none") +
  theme(axis.text = element_text(size = 30, family = "Helvetica"),
        legend.text = element_text(size = 30, family = "Helvetica"),
        legend.title = element_blank(),
        legend.key.size = unit(1.4,"cm")) 
dev.off()

##Heatmap for JGN
JGN <- c("JGN1",
         "JGN2",
         "JGN3",
         "JGN4",
         "JGN5")

neurons <- subset(data, idents = JGN)
DefaultAssay(neurons) = "originalexp"
levels(Idents(neurons))

all.genes <- rownames(neurons)
neurons <- ScaleData(neurons, features = all.genes)

#Average expression for each identity class
JGN <- AverageExpression(neurons,assays = "originalexp",group.by = "ident",return.seurat = T)
#Pull the list of colors used in DimPlot
jugular_col <- c("#00B3F1","#00B0F6","#00ACFB","#00A8FF","#06A4FF")
PalettePlot(palette = jugular_col)

markers <- c("Mrgprd","Gfra3","Trappc3l","Wfdc2","Trpm8")

pdf("figure_heatmap_JGN.pdf", width = 5, height = 7)
DoHeatmap(JGN, features = genes_to_ens(markers), size = 9, raster = F, 
          angle = 90, draw.lines = 0, group.colors = jugular_col) +
  scale_y_discrete(labels = rev(markers)) + 
  scale_fill_viridis()+
  guides(color= "none") +
  theme(axis.text = element_text(size = 25, family = "Helvetica"),
        legend.text = element_text(size = 25, family = "Helvetica"),
        legend.title = element_blank(),
        legend.key.size = unit(1.4,"cm")) 
dev.off()

##Heatmap for non-neuronal cells
NNC <- c("EC1",
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
         "SGC8")

neurons <- subset(data, idents = NNC)
DefaultAssay(neurons) = "originalexp"
levels(Idents(neurons))
all.genes <- rownames(neurons)
neurons <- ScaleData(neurons, features = all.genes)

#Average expression for each identity class
NNC <- AverageExpression(neurons,assays = "originalexp",group.by = "ident",return.seurat = T)
#Pull the list of colors used in DimPlot
nnc_col <- c("#FD6F87","#FC717F","#FA7476","#F8766D","#F47A5F","#F17E4F",
             "#ED813C","#DE8C00","#D89000","#A3A500","#99A800","#8EAB00",
             "#82AD00","#75AF00","#65B200","#53B400","#39B600","#00B70C",
             "#00B931","#00BB45","#00BC56","#00BD64","#00BE71","#00BF7D",
             "#00C089","#00C094")
PalettePlot(palette = nnc_col)

markers <- c("Cldn5",
             "Selp",
             "Podxl",
             "Lum",
             "Acta2",
             "Igfbp6",
             "Bnc2",
             "C1qc",
             "Cd52",
             "Kcna1",
             "Sostdc1",
             "Tyrobp",
             "Ncmap",
             "Prx",
             "Slc36a2",
             "Mt2",
             "Smoc2",
             "Ifitm1",
             "Mmd2",
             "Cxcl10",
             "Bcan",
             "Atf3",
             "Pou3f1",
             "Ifit3",
             "Sdc4",
             "Ttyh1")

pdf("figure_heatmap_NNC.pdf", width = 13, height = 13)
DoHeatmap(NNC, features = genes_to_ens(markers), size = 10, raster = F, 
          angle = 90, draw.lines = 0, group.colors = nnc_col) +
  scale_y_discrete(labels = rev(markers)) + 
  scale_fill_viridis()+
  guides(color= "none") +
  theme(axis.text = element_text(size = 30, family = "Helvetica"),
        legend.text = element_text(size = 30, family = "Helvetica"),
        legend.title = element_blank(),
        legend.key.size = unit(1.4,"cm")) 
dev.off()

#Dotplot
pdf("figure_NGN_dotplot.pdf", width = 11, height = 10)
DotPlot(neurons, features = genes_to_ens(rev(markers)), cols = c(rgb(0.9,0.9,0.9,0.3), "darkred")) + 
  scale_size(range = c(4,10)) +
  scale_x_discrete(labels = rev(markers)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(size = 27, family = "Helvetica"),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 27, family = "Helvetica"),
        legend.spacing = unit(1.0, 'cm'),
        legend.key.size = unit(1.4,"cm")) 
dev.off()

pdf("figure_JGN_dotplot.pdf", width = 11, height = 10)
DotPlot(neurons, features = genes_to_ens(rev(markers)), cols = c(rgb(0.9,0.9,0.9,0.3), "darkred")) + 
  scale_size(range = c(4,10)) +
  scale_x_discrete(labels = rev(markers)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(size = 27, family = "Helvetica"),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 27, family = "Helvetica"),
        legend.spacing = unit(1.0, 'cm'),
        legend.key.size = unit(1.4,"cm")) 
dev.off()

pdf("figure_NNC_dotplot.pdf", width = 11, height = 10)
DotPlot(neurons, features = genes_to_ens(rev(markers)), cols = c(rgb(0.9,0.9,0.9,0.3), "darkred")) + 
  scale_size(range = c(4,10)) +
  scale_x_discrete(labels = rev(markers)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text = element_text(size = 27, family = "Helvetica"),
        axis.title = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 27, family = "Helvetica"),
        legend.spacing = unit(1.0, 'cm'),
        legend.key.size = unit(1.4,"cm")) 
dev.off()

#RCTD
data.st <- readRDS("../../2023_Sep_nodose_spatial/nodose_st_RCTD_2401.RDS")
DefaultAssay(data.st) <- "RNA"
data.st <- SetIdent(data.st, value = "cell_types") 
Idents(data.st) <- factor(Idents(data.st), levels = celltypes)
names(celltype_col)<-levels(data.st)

pdf("figure_RCTD.pdf", width = 3, height = 3)
SpatialDimPlot(data.st, images = "F1R1", stroke = NA, pt.size.factor = 2, cols = celltype_col) + NoLegend()
SpatialDimPlot(data.st, images = "F1R1", stroke = NA, pt.size.factor = 1, cols = celltype_col) + NoLegend()

SpatialDimPlot(data.st, images = "F2R2", stroke = NA, pt.size.factor = 2, cols = celltype_col) + NoLegend()
SpatialDimPlot(data.st, images = "F2R2", stroke = NA, pt.size.factor = 1, cols = celltype_col) + NoLegend()

SpatialDimPlot(data.st, images = "F3R3", stroke = NA, pt.size.factor = 2, cols = celltype_col) + NoLegend()
SpatialDimPlot(data.st, images = "F3R3", stroke = NA, pt.size.factor = 1, cols = celltype_col) + NoLegend()
dev.off()

pdf("figure_st_Phox2b.pdf", width = 3, height = 3)
SpatialFeaturePlot(data.st, "Phox2b", images = "F3R3", stroke = NA, alpha = c(0.1,1), pt.size.factor = 2) +
  NoLegend() +
  scale_fill_gradient(high = "violetred4", low = "grey40", limits=c(0,6))
dev.off()

#Fedfasted
#data$Nutr <- ifelse(data$Nutr.cond == "adlib", "Adlib", "Fast")
pdf("figure_fedfasted.pdf", width = 4, height = 4)
DimPlot(data, group.by = "Nutr", raster = F, order = T,
        cols = c('Adlib'='#117733', 'Fast'='#AA4499'), na.value = rgb(.9,.9,.9,.1)) + 
  ggtitle('') + NoAxes() +
  theme(legend.text = element_text(size = 10, family = "Helvetica")) + coord_fixed()
dev.off()

df <- data.frame(table(data@meta.data$Nutr.cond))
df$Var1 <- c("Adlib", "Fast")

pdf("Figure_fedfasted_barplot.pdf", height = 4, width = 3)
ggplot(df, aes(x = Var1, y = Freq)) +
  geom_bar(stat="identity", fill = c('#117733', '#AA4499'), width = 0.5) + theme_classic() +
  xlab(NULL) +
  ylab("Single Cell / Nucleus") +
  theme(axis.text = element_text(size = 25, color = "black", family = "Helvetica"),
        axis.title = element_text(size = 25, vjust = -4, family = "Helvetica")) 
dev.off()

#Leftright
#data$Po <- ifelse(data$Position == "left", "Left", "Right")
pdf("figure_leftright.pdf", width = 4, height = 4)
DimPlot(data, group.by = "Po", raster = F, order = T, shuffle = T,
        cols = c('Left'='yellowgreen','Right'='plum3'), na.value = rgb(.9,.9,.9,.1)) + 
  ggtitle('') + NoAxes() +
  theme(legend.text = element_text(size = 10, family = "Helvetica")) + coord_fixed()
dev.off()

df <- data.frame(table(data@meta.data$Position))
df$Var1 <- c("Left","Right")

pdf("figure_leftright_barplot.pdf", width = 3, height = 4)
ggplot(df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = c("yellowgreen","plum3"), width = 0.5) + theme_classic() +
  xlab(NULL) +
  ylab("Single Cell / Nucleus") +
  theme(axis.text = element_text(size = 25, colour = "black", family = "Helvetica"),
        axis.title = element_text(size = 25, vjust = -4, family = "Helvetica")) 
dev.off()

#Supplementary data
palette(brewer.pal(7, "Dark2"))

pdf("figure_s_data.pdf", width = 3, height = 3)
FeaturePlot(data, "nCount_originalexp", raster = F, order = T, cols = c(rgb(0.9,0.9,0.9,0.1),palette()[2])) + 
  NoAxes() + ggtitle(NULL) +
  theme(legend.text = element_text(size = 10, family = "Helvetica"),
        legend.spacing = unit(1.0, 'cm')) +
  coord_fixed()

FeaturePlot(data, "nFeature_originalexp", raster = F, cols = c(rgb(0.9,0.9,0.9,0.1),palette()[2])) + 
  NoAxes() + ggtitle(NULL) +
  theme(legend.text = element_text(size = 10, family = "Helvetica"),
        legend.spacing = unit(1.0, 'cm')) +
  coord_fixed()
dev.off()

#Dimplot with different studies
x1 <- DimPlot(data, group.by = "Dataset", raster = FALSE, shuffle = T) + 
  ggtitle("") + NoAxes() +
  scale_color_manual(labels = c("Bai et al., 2019, scSeq",
                                "Buchanan et al., 2022, scSeq",
                                "Cheng et al., 2025, nucSeq",
                                "Kupari et al., 2019, scSeq",
                                "Zhao et al., 2022, scSeq"),
                     values = c("#A50026","#364B9A","#FDB366","#98CAE1","#EAECCC")) + coord_fixed()

pdf("figure_datasets.pdf", width = 4, height = 4)
x1 + NoLegend()
x1 + theme(legend.text = element_text(size = 10, family = "Helvetica"),
        legend.spacing = unit(1.0, 'cm'))
dev.off()

#Feature plots of cell markers
cellmarkers <- c("Emcn","Ebf2","Ptprc","Mpz","Fabp7","Prdm12")
df <- as.data.frame(cellmarkers)
cellmarker_col <- c("#FD6F87","#F17E4F","#D89000","#65B200","#00BE71","#00ACFB")

#palette(brewer.pal(6, "Dark2"))

for (i in 1:length(df$cellmarkers)) {
  scales <- max(summary(GetAssayData(data,slot = "data")[genes_to_ens(df$cellmarkers[i]),]))
  df$scales[i] <- scales
}

pdf("figure_cellmarkers.pdf", width = 3, height = 3)
for (i in 1:length(df$cellmarkers)) {
  p <- FeaturePlot(data, genes_to_ens(df$cellmarkers[i]), raster = F) + 
    ggtitle(NULL) + NoAxes() +
    scale_colour_gradient(high = cellmarker_col[i], low = rgb(0.9,0.9,0.9,0.1), limits=c(0,df$scales[i])) +
    theme(legend.text = element_text(size = 10, family = "Helvetica"),
          legend.spacing = unit(1.0, 'cm')) + 
    coord_fixed()
  plot(p)
}
dev.off()

pdf("figure_st_cellmarkers.pdf", width = 3, height = 3)
for (i in 1:length(df$cellmarkers)) {
  p <- SpatialFeaturePlot(data.st, df$cellmarkers[i], images = "F2R2", stroke = NA, alpha = c(0.1,1), pt.size.factor = 4) +
    NoLegend() +
    scale_fill_gradient(high = cellmarker_col[i], low = "grey60", limits=c(0,df$scales[i]))
  plot(p)
}
dev.off()
























