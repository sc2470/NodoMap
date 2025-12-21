library(Seurat)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(RColorBrewer)
library(viridis)


# load data ---------------------------------------------------------------
data <- readRDS("diet_nodomap_cleaned.RDS")
DefaultAssay(data) = "originalexp"
View(data@meta.data)

data <- SetIdent(data, value = "cluster_number")
levels(Idents(data))

mouse_gene_names<-read.delim("mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}


# Figure 1 ----------------------------------------------------------------
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
Idents(data) <- factor(Idents(data), levels = celltypes)

legend <- c("EC1.Cldn5/Ly6c1/Ly6a",
            "EC2.Selp/Fabp4/Vwf",
            "EC3.Podxl/Depp1/Ptprb",
            "FB1.Lum/Ccl11/Smoc2",
            "FB2.Acta2/Myh11/Ndufa4l2",
            "FB3.Igfbp6/Thbs4/Islr",
            "FB4.Hbb-bs/Bnc2/Ebf2",
            "HC1.C1qa/C1qc/C1qb",
            "HC2.Cd52/Lsp1/Rac2",
            "GC1.Kcna1/Fxyd3/Mbp",
            "GC2.Sostdc1/Plekhb1/Tsc22d4",
            "GC3.Tyrobp/Lyz2/Fcer1g",
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
            "JGN5.Trpm8/Foxp2/Cdh8",
            "NGN1.C1ql3/Adcyap1/St18",
            "NGN2.Htr3b/Htr3a/Gpr65",
            "NGN3.Cysltr2/Chrnb3/Kcnip4",
            "NGN4.Uts2b/Cckar/Vip",
            "NGN5.Kcng1/Vmn1r85/Trpv1",
            "NGN6.Lypd1/Ddc/Tafa2",
            "NGN7.Rbp4/Gda/Kcnk9",
            "NGN8.Trpa1/Efcab6/Nos1",
            "NGN9.Miat/Snhg11/Meg3",
            "NGN10.Olfm3/Tmem233/P2ry1",
            "NGN11.Sprr1a/Ecel1/Ucn",
            "NGN12.Slc18a3/Pappa2/Runx3",
            "NGN13.Lox/Pkib/F2r",
            "NGN14.Slc17a7/Hapln4/Lgi3",
            "NGN15.Bmp3/Gal/Pcdh9",
            "NGN16.Gabra1/Chodl/Gabrb2",
            "NGN17.Thsd7b/Cacng5/Gata3",
            "NGN18.Chrm2/Col24a1/Lrp1b",
            "NGN19.Slc6a2/Npy/Hand2",
            "NGN20.Glp1r/Amigo2/Npy2r",
            "NGN21.Olfr78/Runx3/Avpr1a")

names(legend)<-levels(data)
data<-RenameIdents(data,legend)
data <- AddMetaData(data, metadata = Idents(data), col.name = "cluster_markers")


View(data@meta.data)
#saveRDS(data, "diet_nodomap_Wilcox_2510.rds")

#Set colour scale for different celltypes
celltype_col <- c("#FD6F87","#FC717F","#FA7476","#F8766D","#F47A5F","#F17E4F",
                  "#ED813C","#DE8C00","#D89000","#A3A500","#99A800","#8EAB00",
                  "#82AD00","#75AF00","#65B200","#53B400","#39B600","#00B70C",
                  "#00B931","#00BB45","#00BC56","#00BD64","#00BE71","#00BF7D",
                  "#00C089","#00C094","#00B3F1","#00B0F6","#00ACFB","#00A8FF",
                  "#06A4FF","#8893FF","#9590FF","#A18CFF","#AC88FF","#B584FF",
                  "#BF80FF","#C77CFF","#CF78FF","#D675FE","#DC71FA","#E26EF7",
                  "#E76BF3","#EC69EF","#F066EA","#F464E5","#F763E0","#FA62DB",
                  "#FC61D5","#FE61CF","#FF61C9","#FF61C2")

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
         "NGN21")

Idents(data) <- "cluster_number"
neurons <- subset(data, idents = "NGN")
DefaultAssay(neurons) = "originalexp"
Idents(neurons) <- factor(Idents(neurons), levels = NGN)
levels(Idents(neurons))
all.genes <- rownames(neurons)
neurons <- ScaleData(neurons, features = all.genes)

#Average expression for each identity class
NGN <- AverageExpression(neurons,assays = "originalexp",group.by = "ident",return.seurat = T)
#Pull the list of colors used in DimPlot
nodose_col <- c("#8893FF","#9590FF","#A18CFF","#AC88FF","#B584FF","#BF80FF",
                "#C77CFF","#CF78FF","#D675FE","#DC71FA","#E26EF7","#E76BF3",
                "#EC69EF","#F066EA","#F464E5","#F763E0","#FA62DB","#FC61D5",
                "#FE61CF","#FF61C9","#FF61C2")

markers <- c("C1ql3",
             "Htr3b",
             "Cysltr2",
             "Uts2b",
             "Kcng1",
             "Lypd1",
             "Rbp4",
             "Trpa1",
             "Miat",
             "Olfm3",
             "Sprr1a",
             "Slc18a3",
             "Lox",
             "Slc17a7",
             "Bmp3",
             "Gabra1",
             "Thsd7b",
             "Chrm2",
             "Slc6a2",
             "Glp1r",
             "Olfr78")

pdf("figure_heatmap_NGN.pdf", width = 11, height = 12.5)
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
Idents(neurons) <- factor(Idents(neurons), levels = JGN)
levels(Idents(neurons))
all.genes <- rownames(neurons)
neurons <- ScaleData(neurons, features = all.genes)

#Average expression for each identity class
JGN <- AverageExpression(neurons,assays = "originalexp",group.by = "ident",return.seurat = T)
#Pull the list of colors used in DimPlot
jugular_col <- c("#00B3F1","#00B0F6","#00ACFB","#00A8FF","#06A4FF")

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
Idents(neurons) <- factor(Idents(neurons), levels = NNC)
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

markers <- c("Cldn5",
             "Selp",
             "Podxl",
             "Lum",
             "Acta2",
             "Igfbp6",
             "Hbb-bs",
             "C1qa",
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
             "G0s2",
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


#RCTD
data.st <- readRDS("nodose_st_RCTD_paramrelax_251025.RDS")
DefaultAssay(data.st) <- "RNA"
View(data.st@meta.data)
data.st <- SetIdent(data.st, value = "first_type") 
Idents(data.st) <- factor(Idents(data.st), levels = celltypes)
levels(Idents(data.st))
names(celltype_col)<-levels(data.st)

x<-data.st@meta.data[,c(17,18)]
levels(x$first_type)
x$first_type[x$spot_class != 'singlet'] <- NA
colnames(x)<-c('spot_class', 'singlet_cells')
data.st <- AddMetaData(data.st, metadata = x$singlet_cells, col.name = "singlet_cells")

data.st@images$F1R1@coordinates[,1] <- data.st@images$F1R1@coordinates[,1]*-1
data.st@images$F1R1@coordinates[,2] <- data.st@images$F1R1@coordinates[,2]*-1

pdf("RCTD_Dimplot.pdf", width = 9.6, height = 9.6)
SpatialDimPlot(data.st, group.by = "singlet_cells", images = "F1R1", stroke = NA, 
               pt.size.factor = 1, image.alpha = 0) + NoLegend() + 
  scale_fill_manual(values = celltype_col, na.value="#d6d6d6")

SpatialDimPlot(data.st, group.by = "singlet_cells", images = "F2R2", stroke = NA, 
               pt.size.factor = 1, image.alpha = 0) + NoLegend() + 
  scale_fill_manual(values = celltype_col, na.value="#d6d6d6")

SpatialDimPlot(data.st, group.by = "singlet_cells", images = "F3R3", stroke = NA, 
               pt.size.factor = 1, image.alpha = 0) + NoLegend() + 
  scale_fill_manual(values = celltype_col, na.value="#d6d6d6")
dev.off()


#Feature plots of cell markers
cellmarkers <- c("Emcn","Ebf2","Ptprc","Mpz","Fabp7","Prdm12","Phox2b")
df <- as.data.frame(cellmarkers)
cellmarker_col <- c("#FD6F87","#F17E4F","#D89000","#65B200","#00BE71","#00ACFB","violetred4")

for (i in 1:length(df$cellmarkers)) {
  scales <- max(summary(GetAssayData(data,slot = "data")[genes_to_ens(df$cellmarkers[i]),]))
  df$scales[i] <- scales
}

pdf("figure_cellmarkers.pdf", width = 4, height = 4)
for (i in 1:length(df$cellmarkers)) {
  p <- FeaturePlot(data, genes_to_ens(df$cellmarkers[i]), raster = F, pt.size = 0.1) + 
    ggtitle(NULL) + NoAxes() +
    scale_colour_gradient(high = cellmarker_col[i], low = rgb(0.95,0.95,0.95), limits=c(0,df$scales[i])) +
    theme(legend.text = element_text(size = 10, family = "Helvetica"),
          legend.spacing = unit(1.0, 'cm')) + 
    coord_fixed()
  plot(p)
}
dev.off()




# Figure 2 ----------------------------------------------------------------
#Fedfasted
pdf("figure_fedfasted.pdf", width = 4, height = 4)
DimPlot(data, group.by = "Nutr.cond", raster = F, order = T,
        cols = c('adlib'='#117733', 'fast'='#AA4499'), na.value = rgb(.95,.95,.95)) + 
  ggtitle('') + NoAxes() +
  theme(legend.text = element_text(size = 10, family = "Helvetica")) + coord_fixed()
dev.off()

df <- data.frame(table(data@meta.data$Nutr.cond))
df$Var1 <- c("Fed", "Fast")
df$Var1 <- factor(df$Var1, levels = c("Fed", "Fast"))

pdf("figure_fedfasted_barplot.pdf", height = 4, width = 3)
ggplot(df, aes(x = Var1, y = Freq)) +
  geom_bar(stat="identity", fill = c('#117733', '#AA4499'), width = 0.5) + theme_classic() +
  xlab(NULL) +
  ylab("Single Cell / Nucleus") +
  theme(axis.text = element_text(size = 25, color = "black", family = "Helvetica"),
        axis.title = element_text(size = 25, vjust = -4, family = "Helvetica")) 
dev.off()


# Figure 3 ----------------------------------------------------------------
#Leftright
pdf("figure_leftright.pdf", width = 4, height = 4)
DimPlot(data, group.by = "Position", raster = F, order = T, shuffle = T,
        cols = c('left'='yellowgreen','right'='plum3'), na.value = rgb(.95,.95,.95)) + 
  ggtitle('') + NoAxes() + 
  theme(legend.text = element_text(size = 10, family = "Helvetica")) + coord_fixed()
dev.off()

df <- data.frame(table(data@meta.data$Position))
df$Var1 <- c("Left","Right")

pdf("figure_leftright_barplot.pdf", width = 4, height = 4)
ggplot(df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = c("yellowgreen","plum3"), width = 0.5) + theme_classic() +
  xlab(NULL) +
  ylab("Single Cell / Nucleus") +
  theme(axis.text = element_text(size = 25, colour = "black", family = "Helvetica"),
        axis.title = element_text(size = 25, vjust = -4, family = "Helvetica")) 
dev.off()



# SFigure 1 ---------------------------------------------------------------
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


# SFigure 2 ---------------------------------------------------------------
# Phylogenetic tree -------------------------------------------------------
library(ape)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 10000)

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
         "NGN21")

NG <- subset(data, idents = NGN)
View(NG@meta.data)
levels(Idents(NG))
NG <- SetIdent(NG, value = "cluster_markers")

NG <- BuildClusterTree(NG, reduction = "umap", reorder = T)
Tool(NG, slot = 'BuildClusterTree')
pdf("NGN_Phylogenetic_Tree.pdf", width = 7, height = 5)
plot(NG@tools$BuildClusterTree, type = "phylogram", font = 4, cex = 0.7)
dev.off()

JGN <- c("JGN1",
         "JGN2",
         "JGN3",
         "JGN4",
         "JGN5")
JG <- subset(data, idents = JGN)
JG <- SetIdent(JG, value = "cluster_markers")
JG <- BuildClusterTree(JG, reduction = "umap", reorder = T)
Tool(JG, slot = 'BuildClusterTree')
pdf("JGN_Phylogenetic_Tree.pdf", width = 7, height = 5)
plot(JG@tools$BuildClusterTree, type = "phylogram", font = 4, cex = 0.7)
dev.off()




