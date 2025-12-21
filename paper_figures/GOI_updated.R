library(Seurat)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(RColorBrewer)
library(viridis)

mouse_gene_names<-read.delim("mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

data <- readRDS("../data_objects/diet_nodomap_metadata_updated_2510.rds")
DefaultAssay(data) = "originalexp"
data <- SetIdent(data, value = "cluster_number")
levels(Idents(data))


# FedFast -----------------------------------------------------------------
#Subset fed and fasted data 
inhouse <- subset(data, subset = Dataset == 'inhouse')
inhouse <- SetIdent(inhouse, value = "seurat_clusters")
View(inhouse@meta.data)

#Load the marker list
fedfastmarkers <- list.files(path = ".", pattern = "^fedfasted\\d+.txt")
table(paste0(inhouse$seurat_clusters, inhouse$Nutr.cond))

##some groups have fewer than 3 cells
#7, 8, 18, 19, 20, 30, 31, 33, 34, 39, 42, 45, 46, 48, 49, 50, 51, 52
clustermarkers <- c("EC1","SGC1","SGC2","MGC1","SGC3","NGN1","HC1","FB1","SGC6",
                    "MGC2","FB2","NGN2","GC1","NGN3","MGC3","NGN4","NGN5","SGC8",
                    "JGN1","NGN6","GC2","JGN2","FB3","NGN7","NGN8","NGN10","EC3",
                    "NGN13","NGN14","JGN3","NGN15","JGN4","NGN17","FB4","JGN5")

fedfastmarkers <- do.call(rbind, Map("cbind", lapply(fedfastmarkers, read.delim), cluster=fedfastmarkers))
fedfastmarkers$cluster <- gsub("\\D","",fedfastmarkers$cluster)
colnames(fedfastmarkers)[1] <- "gene"
#order the cluster number
fedfastmarkers$cluster <- as.numeric(fedfastmarkers$cluster)
fedfastmarkers <- fedfastmarkers[order(fedfastmarkers$cluster), ]
unique(fedfastmarkers$cluster)
length(unique(fedfastmarkers$cluster))
#map cell types to cluster numbers
unique_clusters <- sort(unique(fedfastmarkers$cluster))
cluster_to_marker <- setNames(clustermarkers, unique_clusters) #create a named vector
fedfastmarkers$clustermarker <- cluster_to_marker[as.character(fedfastmarkers$cluster)]
unique(paste0(fedfastmarkers$cluster,fedfastmarkers$clustermarker))

fedfastmarkers <- subset(fedfastmarkers, subset = p_val < 0.05)
fedfastmarkers$prop <- fedfastmarkers$avg_log2FC*(fedfastmarkers$pct.1 + 0.01)/(fedfastmarkers$pct.2 + 0.01)

##Subset certain NGN clusters and view their fed&fasted markers
##NGN2.Htr3b/Htr3a/Gpr65
NGN2 <- subset(inhouse, subset = cluster_number == 'NGN2')
DefaultAssay(NGN2) <- "originalexp"
NGN2 <- SetIdent(NGN2, value = "Nutr.cond")
VlnPlot(NGN2, genes_to_ens("Aplp1"))
VlnPlot(NGN2, genes_to_ens("Xkr6"))

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


##NGN10.Olfm3/Tmem233/P2ry1
NGN10 <- subset(inhouse, subset = cluster_number == 'NGN10')
DefaultAssay(NGN10) <- "originalexp"
NGN10 <- SetIdent(NGN10, value = "Nutr.cond")
VlnPlot(NGN10, genes_to_ens("Grm8"))
VlnPlot(NGN10, genes_to_ens("Bbx"))

pdf("figure_fedfast_NGN10_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN10, features = genes_to_ens(c("Grm8","Bbx")))
p[[1]]$labels$y <- "Grm8"
p[[2]]$labels$y <- "Bbx"
levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()


##NGN5.Kcng1/Vmn1r85/Trpv1
NGN5 <- subset(inhouse, subset = cluster_number == 'NGN5')
DefaultAssay(NGN5) <- "originalexp"
NGN5 <- SetIdent(NGN5, value = "Nutr.cond")
VlnPlot(NGN5, genes_to_ens("Prkca"))
VlnPlot(NGN5, genes_to_ens("Rpl6"))

pdf("figure_fedfast_NGN5_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN5, features = genes_to_ens(c("Prkca","Rpl6")))
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


##NGN8.Trpa1/Efcab6/Nos1
NGN8 <- subset(inhouse, subset = cluster_number == 'NGN8')
DefaultAssay(NGN8) <- "originalexp"
NGN8 <- SetIdent(NGN8, value = "Nutr.cond")
VlnPlot(NGN8, genes_to_ens("Srsf2"))
VlnPlot(NGN8, genes_to_ens("Spock1"))

pdf("figure_fedfast_NGN8_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN8, features = genes_to_ens(c("Srsf2","Spock1")))
p[[1]]$labels$y <- "Srsf2"
p[[2]]$labels$y <- "Spock1"
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
JGN1 <- subset(inhouse, subset = cluster_number == 'JGN1')
DefaultAssay(JGN1) <- "originalexp"
JGN1 <- SetIdent(JGN1, value = "Nutr.cond")
VlnPlot(JGN1, genes_to_ens("Lingo2"))
VlnPlot(JGN1, genes_to_ens("Tubb3"))

pdf("figure_fedfast_JGN1_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(JGN1, features = genes_to_ens(c("Lingo2","Tubb3")))
p[[1]]$labels$y <- "Lingo2"
p[[2]]$labels$y <- "Tubb3"
levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()


##JGN5.Trpm8/Foxp2/Cdh8
JGN5 <- subset(inhouse, subset = cluster_number == 'JGN5')
DefaultAssay(JGN5) <- "originalexp"
JGN5 <- SetIdent(JGN5, value = "Nutr.cond")
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
MGC1 <- subset(inhouse, subset = cluster_number == 'MGC1')
DefaultAssay(MGC1) <- "originalexp"
MGC1 <- SetIdent(MGC1, value = "Nutr.cond")
VlnPlot(MGC1, genes_to_ens("Tenm2"))
VlnPlot(MGC1, genes_to_ens("Ssbp2"))

pdf("figure_fedfast_MGC1_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(MGC1, features = genes_to_ens(c("Tenm2","Ssbp2")))
p[[1]]$labels$y <- "Tenm2"
p[[2]]$labels$y <- "Ssbp2"
levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()

##SGC1.Mmd2/Ednrb/Fabp7
SGC1 <- subset(inhouse, subset = cluster_number == 'SGC1')
DefaultAssay(SGC1) <- "originalexp"
SGC1 <- SetIdent(SGC1, value = "Nutr.cond")
VlnPlot(SGC1, genes_to_ens("Adgrg6"))
VlnPlot(SGC1, genes_to_ens("Hmgcs1"))

pdf("figure_fedfast_SGC1_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(SGC1, features = genes_to_ens(c("Adgrg6","Hmgcs1")))
p[[1]]$labels$y <- "Adgrg6"
p[[2]]$labels$y <- "Hmgcs1"
levels(p[[1]]$data$ident) <- c("Adlib", "Fast")
levels(p[[2]]$data$ident) <- c("Adlib", "Fast")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()







# RightLeft ---------------------------------------------------------------
position <- subset(data, subset = Position != "unassigned")
position <- SetIdent(position, value = "seurat_clusters")
View(position@meta.data)

#Load the marker list
leftrightmarkers <- list.files(path = ".", pattern = "^leftright\\d+.txt")
table(paste0(position$seurat_clusters, position$Position))

clustermarker <- c("EC1","SGC1","SGC2","MGC1","SGC3","NGN1","HC1","SGC4","SGC5",
                   "FB1","SGC6","MGC2","FB2","NGN2","GC1","NGN3","MGC3","NGN4",
                   "SGC7","MGC4","NGN5","SGC8","JGN1","NGN6","GC2","JGN2","FB3",
                   "NGN7","NGN8","EC2","NGN9","NGN10","NGN11","NGN12","EC3","NGN13",
                   "NGN14","JGN3","HC2","NGN15","JGN4","NGN16","NGN17","FB4","GC3",
                   "NGN18","JGN5","MGC5","MGC6","NGN19","NGN20","NGN21")

leftrightmarkers <- do.call(rbind, Map("cbind", lapply(leftrightmarkers, read.delim), cluster=leftrightmarkers))
leftrightmarkers$cluster <- gsub("\\D","",leftrightmarkers$cluster)
colnames(leftrightmarkers)[1] <- "gene"
#order the cluster number
leftrightmarkers$cluster <- as.numeric(leftrightmarkers$cluster)
leftrightmarkers <- leftrightmarkers[order(leftrightmarkers$cluster),]
unique(leftrightmarkers$cluster)
length(unique(leftrightmarkers$cluster))
#map cell types to cluster numbers
unique_clusters <- sort(unique(leftrightmarkers$cluster))
cluster_to_marker <- setNames(clustermarker, unique_clusters) #create a named vector
leftrightmarkers$clustermarker <- cluster_to_marker[as.character(leftrightmarkers$cluster)]
unique(paste0(leftrightmarkers$cluster,leftrightmarkers$clustermarker))

leftrightmarkers <- subset(leftrightmarkers, subset = p_val < 0.05)
leftrightmarkers$prop <- (leftrightmarkers$avg_log2FC*(leftrightmarkers$pct.1+0.01))/(leftrightmarkers$pct.2+0.01)

##Subset certain NGN clusters and view their right&left markers
##NGN3.Cysltr2/Chrnb3/Kcnip4
NGN3 <- subset(position, subset = cluster_number == 'NGN3')
DefaultAssay(NGN3) <- "originalexp"
NGN3 <- SetIdent(NGN3, value = "Position")
VlnPlot(NGN3, genes_to_ens("Cntnap2"))
VlnPlot(NGN3, genes_to_ens("Prkg1"))

pdf("figure_leftright_NGN3_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN3, features = genes_to_ens(c("Cntnap2","Prkg1")))
p[[1]]$labels$y <- "Cntnap2"
p[[2]]$labels$y <- "Prkg1"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()

##NGN9.Miat/Snhg11/Meg3
NGN9 <- subset(position, subset = cluster_number == 'NGN9')
DefaultAssay(NGN9) <- "originalexp"
NGN9 <- SetIdent(NGN9, value = "Position")
VlnPlot(NGN9, genes_to_ens("Ctnna2"))
VlnPlot(NGN9, genes_to_ens("Dbi"))

pdf("figure_leftright_NGN9_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN9, features = genes_to_ens(c("Ctnna2","Dbi")))
p[[1]]$labels$y <- "Ctnna2"
p[[2]]$labels$y <- "Dbi"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()

##NGN16.Gabra1/Chodl/Gabrb2
NGN16 <- subset(position, subset = cluster_number == 'NGN16')
DefaultAssay(NGN16) <- "originalexp"
NGN16 <- SetIdent(NGN16, value = "Position")
VlnPlot(NGN16, genes_to_ens("Negr1"))
VlnPlot(NGN16, genes_to_ens("Fos"))

pdf("figure_leftright_NGN16_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN16, features = genes_to_ens(c("Negr1","Fos")))
p[[1]]$labels$y <- "Negr1"
p[[2]]$labels$y <- "Fos"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()

##NGN18.Chrm2/Col24a1/Lrp1b
NGN18 <- subset(position, subset = cluster_number == 'NGN18')
DefaultAssay(NGN18) <- "originalexp"
NGN18 <- SetIdent(NGN18, value = "Position")
VlnPlot(NGN18, genes_to_ens("Pcdh9"))
VlnPlot(NGN18, genes_to_ens("Nop10"))

pdf("figure_leftright_NGN18_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(NGN18, features = genes_to_ens(c("Pcdh9","Nop10")))
p[[1]]$labels$y <- "Pcdh9"
p[[2]]$labels$y <- "Nop10"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()

##Subset certain JGN clusters and view their left&right markers
##JGN1.Mrgprd/Tmem45b/Grik1
JGN1 <- subset(position, subset = cluster_number == 'JGN1')
DefaultAssay(JGN1) <- "originalexp"
JGN1 <- SetIdent(JGN1, value = "Position")
VlnPlot(JGN1, genes_to_ens("Htr4"))
VlnPlot(JGN1, genes_to_ens("Actb"))

pdf("figure_leftright_JGN1_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(JGN1, features = genes_to_ens(c("Htr4","Actb")))
p[[1]]$labels$y <- "Htr4"
p[[2]]$labels$y <- "Actb"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()

##JGN3.Trappc3l/Nptx1/Tafa1
JGN3 <- subset(position, subset = cluster_number == 'JGN3')
DefaultAssay(JGN3) <- "originalexp"
JGN3 <- SetIdent(JGN3, value = "Position")
VlnPlot(JGN3, genes_to_ens("Emb"))
VlnPlot(JGN3, genes_to_ens("Plxna2"))

pdf("figure_leftright_JGN3_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(JGN3, features = genes_to_ens(c("Emb","Plxna2")))
p[[1]]$labels$y <- "Emb"
p[[2]]$labels$y <- "Plxna2"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()

##Subset certain NNC clusters and view their left&right markers
##FB3.Igfbp6/Thbs4/Islr
FB3 <- subset(position, subset = cluster_number == 'FB3')
DefaultAssay(FB3) <- "originalexp"
FB3 <- SetIdent(FB3, value = "Position")
VlnPlot(FB3, genes_to_ens("Gpc6"))
VlnPlot(FB3, genes_to_ens("Rps27"))

pdf("figure_leftright_FB3_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(FB3, features = genes_to_ens(c("Gpc6","Rps27")))
p[[1]]$labels$y <- "Gpc6"
p[[2]]$labels$y <- "Rps27"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()

##FB4.Hbb-bs/Bnc2/Ebf2
FB4 <- subset(position, subset = cluster_number == 'FB4')
DefaultAssay(FB4) <- "originalexp"
FB4 <- SetIdent(FB4, value = "Position")
VlnPlot(FB4, genes_to_ens("Hspa8"))
VlnPlot(FB4, genes_to_ens("Pitpnc1"))

pdf("figure_leftright_FB4_GOI.pdf", width = 3, height = 3)
p <- Stacked_VlnPlot(FB4, features = genes_to_ens(c("Hspa8","Pitpnc1")))
p[[1]]$labels$y <- "Hspa8"
p[[2]]$labels$y <- "Pitpnc1"
levels(p[[1]]$data$ident) <- c("Left", "Right")
levels(p[[2]]$data$ident) <- c("Left", "Right")
p[[1]]$theme$text$family <- "Helvetica"
p[[2]]$theme$text$family <- "Helvetica"
p[[1]]$theme$axis.text$size <- 14
p[[2]]$theme$axis.text$size <-14
plot(p)
dev.off()





