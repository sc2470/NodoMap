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

#Subset left and right data 
position <- subset(data, subset = Position != "unassigned")
View(position@meta.data)

#Load the marker list
seurat <- list.files(path = ".", pattern = "^leftright\\d+.txt")
table(paste0(position$integrated_snn_res.1.4, position$Position))

celltype <- c("EC1",
              "SGC1",
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
              "SGC3",
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
              "NGN1",
              "NGN20",
              "NGN21",
              "NGN22",
              "HC1",
              "SGC4",
              "SGC5",
              "FB1")
  
seurat <- do.call(rbind, Map("cbind", lapply(seurat, read.delim), cluster=seurat, type=celltype))
seurat$cluster <- gsub("\\D","",seurat$cluster)
colnames(seurat)[1] <- "gene"
seurat <- subset(seurat, subset = p_val < 0.05)
seurat$markergenes <- (seurat$avg_log2FC*(seurat$pct.1+0.01))/(seurat$pct.2+0.01)

##Subset certain NGN clusters and view their left&right markers
##NGN5.Rtn1/Ngfr/Clu
NGN5 <- subset(position, subset = cell_types == 'NGN5')
DefaultAssay(NGN5) <- "originalexp"
NGN5 <- SetIdent(NGN5, value = "Position")
NGN5_marker <- subset(seurat, subset = type == 'NGN5')
VlnPlot(NGN5, genes_to_ens("Mbp"))
VlnPlot(NGN5, genes_to_ens("Sgcd"))

pdf("figure_leftright_NGN5_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN5, features = genes_to_ens(c("Mbp","Sgcd")))
p[[1]]$labels$y <- "Mbp"
p[[2]]$labels$y <- "Sgcd"
p[[1]]$theme$axis.title.y$face <- "italic"
p[[2]]$theme$axis.title.y$face <- "italic"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14

extract_max <- function(x){
  ymax<-max(x$data[[1]])
  return(ceiling(ymax))
}

ymax <- extract_max(p[[1]])
p[[1]] <- p[[1]] + scale_y_continuous(breaks = c(0,ymax))
ymax <- extract_max(p[[2]])
p[[2]] <- p[[2]] + scale_y_continuous(breaks = c(0,ymax))
plot(p)
dev.off()

##NGN10.Miat/Snhg11/Meg3
NGN10 <- subset(position, subset = cell_types == 'NGN10')
DefaultAssay(NGN10) <- "originalexp"
NGN10 <- SetIdent(NGN10, value = "Position")
NGN10_marker <- subset(seurat, subset = type == 'NGN10')
VlnPlot(NGN10, genes_to_ens("Ctnna2"))
VlnPlot(NGN10, genes_to_ens("Dbi"))

pdf("figure_leftright_NGN10_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN10, features = genes_to_ens(c("Ctnna2","Dbi")))
p[[1]]$labels$y <- "Ctnna2"
p[[2]]$labels$y <- "Dbi"
p[[1]]$theme$axis.title.y$face <- "italic"
p[[2]]$theme$axis.title.y$face <- "italic"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
ymax <- extract_max(p[[1]])
p[[1]] <- p[[1]] + scale_y_continuous(breaks = c(0,ymax))
ymax <- extract_max(p[[2]])
p[[2]] <- p[[2]] + scale_y_continuous(breaks = c(0,ymax))
plot(p)
dev.off()

##NGN17.Gabra1/Chodl/Gabrb2
NGN17 <- subset(position, subset = cell_types == 'NGN17')
DefaultAssay(NGN17) <- "originalexp"
NGN17 <- SetIdent(NGN17, value = "Position")
NGN17_marker <- subset(seurat, subset = type == 'NGN17')
VlnPlot(NGN17, genes_to_ens("Negr1"))
VlnPlot(NGN17, genes_to_ens("Fos"))

pdf("figure_leftright_NGN17_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN17, features = genes_to_ens(c("Negr1","Fos")))
p[[1]]$labels$y <- "Negr1"
p[[2]]$labels$y <- "Fos"
p[[1]]$theme$axis.title.y$face <- "italic"
p[[2]]$theme$axis.title.y$face <- "italic"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
ymax <- extract_max(p[[1]])
p[[1]] <- p[[1]] + scale_y_continuous(breaks = c(0,ymax))
ymax <- extract_max(p[[2]])
p[[2]] <- p[[2]] + scale_y_continuous(breaks = c(0,ymax))
plot(p)
dev.off()

##NGN19.Gabra1/Gabrb2/Chrm2
NGN19 <- subset(position, subset = cell_types == 'NGN19')
DefaultAssay(NGN19) <- "originalexp"
NGN19 <- SetIdent(NGN19, value = "Position")
NGN19_marker <- subset(seurat, subset = type == 'NGN19')
VlnPlot(NGN19, genes_to_ens("Pcdh9"))
VlnPlot(NGN19, genes_to_ens("Nop10"))

pdf("figure_leftright_NGN19_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN19, features = genes_to_ens(c("Pcdh9","Nop10")))
p[[1]]$labels$y <- "Pcdh9"
p[[2]]$labels$y <- "Nop10"
p[[1]]$theme$axis.title.y$face <- "italic"
p[[2]]$theme$axis.title.y$face <- "italic"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
ymax <- extract_max(p[[1]])
p[[1]] <- p[[1]] + scale_y_continuous(breaks = c(0,ymax))
ymax <- extract_max(p[[2]])
p[[2]] <- p[[2]] + scale_y_continuous(breaks = c(0,ymax))
plot(p)
dev.off()

##Subset certain JGN clusters and view their left&right markers
##JGN1.Mrgprd/Tmem45b/Grik1
JGN1 <- subset(position, subset = cell_types == 'JGN1')
DefaultAssay(JGN1) <- "originalexp"
JGN1 <- SetIdent(JGN1, value = "Position")
JGN1_marker <- subset(seurat, subset = type == 'JGN1')
VlnPlot(JGN1, genes_to_ens("Htr4"))
VlnPlot(JGN1, genes_to_ens("Actb"))

pdf("figure_leftright_JGN1_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(JGN1, features = genes_to_ens(c("Htr4","Actb")))
p[[1]]$labels$y <- "Htr4"
p[[2]]$labels$y <- "Actb"
p[[1]]$theme$axis.title.y$face <- "italic"
p[[2]]$theme$axis.title.y$face <- "italic"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
ymax <- extract_max(p[[1]])
p[[1]] <- p[[1]] + scale_y_continuous(breaks = c(0,ymax))
ymax <- extract_max(p[[2]])
p[[2]] <- p[[2]] + scale_y_continuous(breaks = c(0,ymax))
plot(p)
dev.off()

##JGN3.Trappc3l/Nptx1/Tafa1
JGN3 <- subset(position, subset = cell_types == 'JGN3')
DefaultAssay(JGN3) <- "originalexp"
JGN3 <- SetIdent(JGN3, value = "Position")
JGN3_marker <- subset(seurat, subset = type == 'JGN3')
VlnPlot(JGN3, genes_to_ens("Emb"))
VlnPlot(JGN3, genes_to_ens("Plxna2"))

pdf("figure_leftright_JGN3_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(JGN3, features = genes_to_ens(c("Emb","Plxna2")))
p[[1]]$labels$y <- "Emb"
p[[2]]$labels$y <- "Plxna2"
p[[1]]$theme$axis.title.y$face <- "italic"
p[[2]]$theme$axis.title.y$face <- "italic"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
ymax <- extract_max(p[[1]])
p[[1]] <- p[[1]] + scale_y_continuous(breaks = c(0,ymax))
ymax <- extract_max(p[[2]])
p[[2]] <- p[[2]] + scale_y_continuous(breaks = c(0,ymax))
plot(p)
dev.off()

##Subset certain NNC clusters and view their left&right markers
##FB3.Igfbp6/Thbs4/Islr
FB3 <- subset(position, subset = cell_types == 'FB3')
DefaultAssay(FB3) <- "originalexp"
FB3 <- SetIdent(FB3, value = "Position")
FB3_marker <- subset(seurat, subset = type == 'FB3')
VlnPlot(FB3, genes_to_ens("Gpc6"))
VlnPlot(FB3, genes_to_ens("Rps27"))

pdf("figure_leftright_FB3_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(FB3, features = genes_to_ens(c("Gpc6","Rps27")))
p[[1]]$labels$y <- "Gpc6"
p[[2]]$labels$y <- "Rps27"
p[[1]]$theme$axis.title.y$face <- "italic"
p[[2]]$theme$axis.title.y$face <- "italic"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
ymax <- extract_max(p[[1]])
p[[1]] <- p[[1]] + scale_y_continuous(breaks = c(0,ymax))
ymax <- extract_max(p[[2]])
p[[2]] <- p[[2]] + scale_y_continuous(breaks = c(0,ymax))
plot(p)
dev.off()

##FB4.Hbb-bs/Bnc2/Ebf2
FB4 <- subset(position, subset = cell_types == 'FB4')
DefaultAssay(FB4) <- "originalexp"
FB4 <- SetIdent(FB4, value = "Position")
FB4_marker <- subset(seurat, subset = type == 'FB4')
VlnPlot(FB4, genes_to_ens("Hspa8"))
VlnPlot(FB4, genes_to_ens("Pitpnc1"))

pdf("figure_leftright_FB4_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(FB4, features = genes_to_ens(c("Hspa8","Pitpnc1")))
p[[1]]$labels$y <- "Hspa8"
p[[2]]$labels$y <- "Pitpnc1"
p[[1]]$theme$axis.title.y$face <- "italic"
p[[2]]$theme$axis.title.y$face <- "italic"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
ymax <- extract_max(p[[1]])
p[[1]] <- p[[1]] + scale_y_continuous(breaks = c(0,ymax))
ymax <- extract_max(p[[2]])
p[[2]] <- p[[2]] + scale_y_continuous(breaks = c(0,ymax))
plot(p)
dev.off()


