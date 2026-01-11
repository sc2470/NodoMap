library(Seurat)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(stringr)

#Use the diet Seurat object
data <- readRDS("diet_nodomap_metadata_updated_2510.rds")
DefaultAssay(data)
data <- SetIdent(data, value = "cluster_number")
table(Idents(data))
levels(Idents(data))

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


# Search by gene names-----------------------------------------------------
mouse_gene_names<-read.delim("mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol)+   
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

#plots to visualize interested genes
data_position <- subset(data, subset = Position != "unassigned")
View(data_position@meta.data)

#Gut hormone receptors
GHR <- c("Cckar","Cckbr","Cnr1","Cnr2","Calcrl","Gcgr","Ghsr","Glp1r","Glp2r","Gipr",
         "Insr","Lepr","Npy1r","Npy2r","Npy4r","Npy5r","Ntsr1","Ntsr2","Oxtr","Sstr1",
         "Sstr5")

pdf("figure_LR_GHR.pdf", width = 12, height = 10)
p <- Stacked_VlnPlot(data_position, features = genes_to_ens(GHR), split.by = "Position",
                     x_lab_rotate = TRUE, plot_legend = TRUE)
p[[1]]$labels$y <- "Cckar"
p[[2]]$labels$y <- "Cckbr"
p[[3]]$labels$y <- "Cnr1"
p[[4]]$labels$y <- "Cnr2"
p[[5]]$labels$y <- "Calcrl"
p[[6]]$labels$y <- "Gcgr"
p[[7]]$labels$y <- "Ghsr"
p[[8]]$labels$y <- "Glp1r"
p[[9]]$labels$y <- "Glp2r"
p[[10]]$labels$y <- "Gipr"
p[[11]]$labels$y <- "Insr"
p[[12]]$labels$y <- "Lepr"
p[[13]]$labels$y <- "Npy1r"
p[[14]]$labels$y <- "Npy2r"
p[[15]]$labels$y <- "Npy4r"
p[[16]]$labels$y <- "Npy5r"
p[[17]]$labels$y <- "Ntsr1"
p[[18]]$labels$y <- "Ntsr2"
p[[19]]$labels$y <- "Oxtr"
p[[20]]$labels$y <- "Sstr1"
p[[21]]$labels$y <- "Sstr5"

plot(p)
dev.off()


#Nutrient receptors
NTR <- c("Casr","Ffar1","Ffar2","Ffar3","Ffar4","Gpbar1","Gpr35","Gpr84",
         "Gpr119","Gpr142","Hcar1","Hcar2","Sucnr1","Slc5a1","Slc5a2","Slc5a3",
         "Slc5a5","Slc5a6","Slc5a7","Tas1r1","Tas1r2","Tas1r3")

#plots to visualize interested genes

pdf("figure_LR_NTR.pdf", width = 12, height = 10)
p <- Stacked_VlnPlot(data_position, features = genes_to_ens(NTR), split.by = "Position",
                     x_lab_rotate = TRUE, plot_legend = TRUE)

p[[1]]$labels$y <- "Casr"
p[[2]]$labels$y <- "Ffar1"
p[[3]]$labels$y <- "Ffar2"
p[[4]]$labels$y <- "Ffar3"
p[[5]]$labels$y <- "Ffar4"
p[[6]]$labels$y <- "Gpbar1"
p[[7]]$labels$y <- "Gpr35"
p[[8]]$labels$y <- "Gpr84"
p[[9]]$labels$y <- "Gpr119"
p[[10]]$labels$y <- "Gpr142"
p[[11]]$labels$y <- "Hcar1"
p[[12]]$labels$y <- "Hcar2"
p[[13]]$labels$y <- "Sucnr1"
p[[14]]$labels$y <- "Slc5a1"
p[[15]]$labels$y <- "Slc5a2"
p[[16]]$labels$y <- "Slc5a3"
p[[17]]$labels$y <- "Slc5a5"
p[[18]]$labels$y <- "Slc5a6"
p[[19]]$labels$y <- "Slc5a7"
p[[20]]$labels$y <- "Tas1r1"
p[[21]]$labels$y <- "Tas1r2"
p[[22]]$labels$y <- "Tas1r3"

plot(p)
dev.off()


# FedFasted ---------------------------------------------------------------
data_nutrient <- subset(data, subset = Dataset == 'inhouse')
View(data_nutrient@meta.data)
levels(Idents(data_nutrient))

pdf("figure_FF_GHR.pdf", width = 12, height = 10)
p <- Stacked_VlnPlot(data_nutrient, features = genes_to_ens(GHR), split.by = "Nutr.cond",
                     x_lab_rotate = TRUE, plot_legend = TRUE)
p[[1]]$labels$y <- "Cckar"
p[[2]]$labels$y <- "Cckbr"
p[[3]]$labels$y <- "Cnr1"
p[[4]]$labels$y <- "Cnr2"
p[[5]]$labels$y <- "Calcrl"
p[[6]]$labels$y <- "Gcgr"
p[[7]]$labels$y <- "Ghsr"
p[[8]]$labels$y <- "Glp1r"
p[[9]]$labels$y <- "Glp2r"
p[[10]]$labels$y <- "Gipr"
p[[11]]$labels$y <- "Insr"
p[[12]]$labels$y <- "Lepr"
p[[13]]$labels$y <- "Npy1r"
p[[14]]$labels$y <- "Npy2r"
p[[15]]$labels$y <- "Npy4r"
p[[16]]$labels$y <- "Npy5r"
p[[17]]$labels$y <- "Ntsr1"
p[[18]]$labels$y <- "Ntsr2"
p[[19]]$labels$y <- "Oxtr"
p[[20]]$labels$y <- "Sstr1"
p[[21]]$labels$y <- "Sstr5"

plot(p)
dev.off()

pdf("figure_FF_NTR.pdf", width = 12, height = 10)
p <- Stacked_VlnPlot(data_nutrient, features = genes_to_ens(NTR), split.by = "Nutr.cond",
                     x_lab_rotate = TRUE, plot_legend = TRUE)

p[[1]]$labels$y <- "Casr"
p[[2]]$labels$y <- "Ffar1"
p[[3]]$labels$y <- "Ffar2"
p[[4]]$labels$y <- "Ffar3"
p[[5]]$labels$y <- "Ffar4"
p[[6]]$labels$y <- "Gpbar1"
p[[7]]$labels$y <- "Gpr35"
p[[8]]$labels$y <- "Gpr84"
p[[9]]$labels$y <- "Gpr119"
p[[10]]$labels$y <- "Gpr142"
p[[11]]$labels$y <- "Hcar1"
p[[12]]$labels$y <- "Hcar2"
p[[13]]$labels$y <- "Sucnr1"
p[[14]]$labels$y <- "Slc5a1"
p[[15]]$labels$y <- "Slc5a2"
p[[16]]$labels$y <- "Slc5a3"
p[[17]]$labels$y <- "Slc5a5"
p[[18]]$labels$y <- "Slc5a6"
p[[19]]$labels$y <- "Slc5a7"
p[[20]]$labels$y <- "Tas1r1"
p[[21]]$labels$y <- "Tas1r2"
p[[22]]$labels$y <- "Tas1r3"

plot(p)
dev.off()







