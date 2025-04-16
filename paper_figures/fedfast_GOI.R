library(Seurat)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(RColorBrewer)
library(viridis)

mouse_gene_names<-read.delim("../../mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

data <- readRDS("../figures/nodose_sc_integrated_Wicox_2401.rds")
DefaultAssay(data)
table(Idents(data))

#Subset fed and fasted data 
inhouse <- subset(data, subset = Dataset == 'inhouse')
View(inhouse@meta.data)

#Load the marker list
markers <- list.files(path = ".", pattern = "^fedfasted\\d+.txt")
table(paste0(inhouse$integrated_snn_res.1.4, inhouse$Nutr.cond))

##some groups have fewer than 3 cells
##7, 8, 18, 19, 30, 31, 33, 34, 39, 42, 45, 46, 48, 49, 50, 52 
cells <- c("EC1",
           "SGC1",
           "SGC6",
           "MGC2",
           "FB2",
           "NGN2",
           "GC1",
           "NGN3",
           "MGC3",
           "NGN4",
           "SGC2",
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
           "MGC1",
           "NGN11",
           "EC3",
           "NGN14",
           "NGN15",
           "JGN3",
           "SGC3",
           "NGN16",
           "JGN4",
           "NGN18",
           "FB4",
           "JGN5",
           "NGN1",
           "HC1",
           "FB1")

markers <- do.call(rbind, Map("cbind", lapply(markers, read.delim), cluster=markers, type = cells))
markers$cluster <- gsub("\\D","",markers$cluster)
colnames(markers)[1] <- "gene"
markers <- subset(markers, subset = p_val < 0.05)

##Subset certain NGN clusters and view their fed&fasted markers
##NGN2.Htr3b/Htr3a/Hs3st4
NGN2 <- subset(inhouse, subset = cell_types == 'NGN2')
DefaultAssay(NGN2) <- "originalexp"
NGN2 <- SetIdent(NGN2, value = "Nutr.cond")
NGN2_marker <- subset(markers, subset = type == 'NGN2')
VlnPlot(NGN2, genes_to_ens("Xkr6"))
VlnPlot(NGN2, genes_to_ens("Aplp1"))

pdf("figure_fedfast_NGN2_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN2, features = genes_to_ens(c("Xkr6","Aplp1")))
p[[1]]$labels$y <- "Xkr6"
p[[2]]$labels$y <- "Aplp1"

levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")

p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"

p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14

plot(p)
dev.off()

##NGN6.Kcng1/Vmn1r85/Trpv1
NGN6 <- subset(inhouse, subset = cell_types == 'NGN6')
DefaultAssay(NGN6) <- "originalexp"
NGN6 <- SetIdent(NGN6, value = "Nutr.cond")
NGN6_marker <- subset(markers, subset = cluster == "NGN6")
VlnPlot(NGN6, genes_to_ens("Prkca"))
VlnPlot(NGN6, genes_to_ens("Rpl6"))

pdf("figure_fedfast_NGN6_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN6, features = genes_to_ens(c("Prkca","Rpl6")))
p[[1]]$labels$y <- "Prkca"
p[[2]]$labels$y <- "Rpl6"

levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")

p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"

p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14

plot(p)
dev.off()

##NGN9.Trpa1/Nos1/Efcab6
NGN9 <- subset(inhouse, subset = cell_types == 'NGN9')
DefaultAssay(NGN9) <- "originalexp"
NGN9 <- SetIdent(NGN9, value = "Nutr.cond")
NGN9_marker <- subset(markers, subset = cluster == "NGN9")
VlnPlot(NGN9, genes_to_ens("Srsf2"))
VlnPlot(NGN9, genes_to_ens("Rpl3"))

pdf("figure_fedfast_NGN9_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN9, features = genes_to_ens(c("Srsf2","Rpl3")))
p[[1]]$labels$y <- "Srsf2"
p[[2]]$labels$y <- "Rpl3"

levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")

p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"

p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14

plot(p)
dev.off()

##NGN11.Olfm3/Tmem233/P2ry1
NGN11 <- subset(inhouse, subset = cell_types == 'NGN11')
DefaultAssay(NGN11) <- "originalexp"
NGN11 <- SetIdent(NGN11, value = "Nutr.cond")
NGN11_marker <- subset(markers, subset = cluster == "NGN11")
VlnPlot(NGN11, genes_to_ens("Cdh12"))
VlnPlot(NGN11, genes_to_ens("Bbx"))

pdf("figure_fedfast_NGN11_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN11, features = genes_to_ens(c("Cdh12","Bbx")))
p[[1]]$labels$y <- "Cdh12"
p[[2]]$labels$y <- "Bbx"

levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")

p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"

p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14

plot(p)
dev.off()


##Subset certain JGN clusters and view their fed&fasted markers
##JGN1.Mrgprd/Tmem45b/Grik1
JGN1 <- subset(inhouse, subset = cell_types == 'JGN1')
DefaultAssay(JGN1) <- "originalexp"
JGN1 <- SetIdent(JGN1, value = "Nutr.cond")
JGN1_marker <- subset(markers, subset = type == 'JGN1')
VlnPlot(JGN1, genes_to_ens("Exoc6"))
VlnPlot(JGN1, genes_to_ens("Rpl39"))

pdf("figure_fedfast_JGN1_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(JGN1, features = genes_to_ens(c("Exoc6","Rpl39")))
p[[1]]$labels$y <- "Exoc6"
p[[2]]$labels$y <- "Rpl39"

levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")

p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"

p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14

plot(p)
dev.off()


##JGN5.Trpm8/Foxp2/Gsg1l
JGN5 <- subset(inhouse, subset = cell_types == 'JGN5')
DefaultAssay(JGN5) <- "originalexp"
JGN5 <- SetIdent(JGN5, value = "Nutr.cond")
JGN5_marker <- subset(markers, subset = type == 'JGN5')
VlnPlot(JGN5, genes_to_ens("Grik4"))
VlnPlot(JGN5, genes_to_ens("Prkg2"))

pdf("figure_fedfast_JGN5_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(JGN5, features = genes_to_ens(c("Grik4","Prkg2")))
p[[1]]$labels$y <- "Grik4"
p[[2]]$labels$y <- "Prkg2"

levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")

p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"

p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14

plot(p)
dev.off()

##MGC1.Ncmap/Fam178b/Cldn19
MGC1 <- subset(inhouse, subset = cell_types == 'MGC1')
DefaultAssay(MGC1) <- "originalexp"
MGC1 <- SetIdent(MGC1, value = "Nutr.cond")
MGC1_marker <- subset(markers, subset = type == 'MGC1')
VlnPlot(MGC1, genes_to_ens("Plekhg1"))
VlnPlot(MGC1, genes_to_ens("Zbtb16"))

pdf("figure_fedfast_MGC1_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(MGC1, features = genes_to_ens(c("Plekhg1","Zbtb16")))
p[[1]]$labels$y <- "Plekhg1"
p[[2]]$labels$y <- "Zbtb16"

levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")

p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"

p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14

plot(p)
dev.off()

##SGC1.Mmd2/Ednrb/Fabp7
SGC1 <- subset(inhouse, subset = cell_types == 'SGC1')
DefaultAssay(SGC1) <- "originalexp"
SGC1 <- SetIdent(SGC1, value = "Nutr.cond")
SGC1_marker <- subset(markers, subset = type == 'SGC1')
VlnPlot(SGC1, genes_to_ens("Creb5"))
VlnPlot(SGC1, genes_to_ens("Brip1os"))

pdf("figure_fedfast_SGC1_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(SGC1, features = genes_to_ens(c("Creb5","Brip1os")))
p[[1]]$labels$y <- "Creb5"
p[[2]]$labels$y <- "Brip1os"

levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")

p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"

p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14

plot(p)
dev.off()




