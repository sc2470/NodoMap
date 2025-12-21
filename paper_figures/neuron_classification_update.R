library(Seurat)
library(ggplot2)
library(scCustomize)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(stringr)
library(pheatmap)

# load data ---------------------------------------------------------------
data <- readRDS("diet_nodomap_metadata_updated_2510.rds")
DefaultAssay(data) = "originalexp"

data <- SetIdent(data, value = "cluster_number")
levels(Idents(data))

mouse_gene_names<-read.delim("mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

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

# Subset neuron clusters-------------------------------------------------------------------------
neurons <- c("JGN1",
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
data.n <- subset(data, idents = neurons)
levels(Idents(data.n))

# Total number of cells per cluster ---------------------------------------
all_cell <- data.frame(table(data.n@active.ident))

###Sodium channel
table(GetAssayData(data.n, slot = "data")[genes_to_ens("Scn1a"),]>0)
Scn1a <- subset(data.n, subset = ENSMUSG00000064329 >0)
Scn1a_cell <- data.frame(table(Idents(Scn1a)))
all_cell$Scn1a <- Scn1a_cell$Freq[match(all_cell$Var1, Scn1a_cell$Var1)]
all_cell$Scn1a_prop <- all_cell$Scn1a/all_cell$Freq
hist(all_cell$Scn1a_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Scn1a"), label = T, order = T)
true_Scn1a <- subset(all_cell, subset = all_cell$Scn1a_prop>=0.4)
write.table(true_Scn1a$Var1, "Nav1.1_pc50_res1.4.txt", sep = "\t", quote = F)


table(GetAssayData(data.n, slot = "data")[genes_to_ens("Scn10a"),]>0)
Scn10a <- subset(data.n, subset = ENSMUSG00000034533 >0)
Scn10a_cell <- data.frame(table(Idents(Scn10a)))
all_cell$Scn10a <- Scn10a_cell$Freq[match(all_cell$Var1, Scn10a_cell$Var1)]
all_cell$Scn10a_prop <- all_cell$Scn10a/all_cell$Freq
hist(all_cell$Scn10a_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Scn10a"), label = T, order = T, repel = T)
true_Scn10a <- subset(all_cell, subset = all_cell$Scn10a_prop>=0.3)
write.table(true_Scn10a$Var1, "Nav1.8_pc50_res1.4.txt", sep = "\t", quote = F)


# -------------------------------------------------------------------------
#Myelinated vs Lightly/unmyelinated neurons
##Myelinated A-fibre structure (Nefh, Cntn1, Ncam1(low))
##Thinly myelinated --> combined expression of Nefh, Cntn1, Ncam1
table(GetAssayData(data.n, slot = "data")[genes_to_ens("Nefh"),]>0)
Nefh <- subset(data.n, subset = ENSMUSG00000020396 >0)
Nefh_cell <- data.frame(table(Idents(Nefh)))
all_cell$Nefh <- Nefh_cell$Freq[match(all_cell$Var1, Nefh_cell$Var1)]
all_cell$Nefh_prop <- all_cell$Nefh/all_cell$Freq
hist(all_cell$Nefh_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Nefh"), label = T)
true_Nefh <- subset(all_cell, subset = all_cell$Nefh_prop>=0.8)
write.table(true_Nefh$Var1, "Myelinated_Nefh_pc50_res1.4.txt", sep = "\t", quote = F)


table(GetAssayData(data.n, slot = "data")[genes_to_ens("Cntn1"),]>0)
Cntn1 <- subset(data.n, subset = ENSMUSG00000055022 >0)
Cntn1_cell <- data.frame(table(Idents(Cntn1)))
all_cell$Cntn1 <- Cntn1_cell$Freq[match(all_cell$Var1, Cntn1_cell$Var1)]
all_cell$Cntn1_prop <- all_cell$Cntn1/all_cell$Freq
hist(all_cell$Cntn1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Cntn1"), label = T)
true_Cntn1 <- subset(all_cell, subset = all_cell$Cntn1_prop>=0.8)
write.table(true_Cntn1$Var1, "Myelinated_Cntn1_pc50_res1.4.txt", sep = "\t", quote = F)


table(GetAssayData(data.n, slot = "data")[genes_to_ens("Cntnap1"),]>0)
Cntnap1 <- subset(data.n, subset = ENSMUSG00000017167 >0)
Cntnap1_cell <- data.frame(table(Idents(Cntnap1)))
all_cell$Cntnap1 <- Cntnap1_cell$Freq[match(all_cell$Var1, Cntnap1_cell$Var1)]
all_cell$Cntnap1_prop <- all_cell$Cntnap1/all_cell$Freq
hist(all_cell$Cntnap1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Cntnap1"), label = T, order = T)
true_Cntnap1 <- subset(all_cell, subset = all_cell$Cntnap1_prop>=0.5)
write.table(true_Cntnap1$Var1, "Myelinated_Cntnap1_pc50_res1.4.txt", sep = "\t", quote = F)


table(GetAssayData(data.n, slot = "data")[genes_to_ens("Ncam1"),]>0)
Ncam1 <- subset(data.n, subset = ENSMUSG00000039542 >0)
Ncam1_cell <- data.frame(table(Idents(Ncam1)))
all_cell$Ncam1 <- Ncam1_cell$Freq[match(all_cell$Var1, Ncam1_cell$Var1)]
all_cell$Ncam1_prop <- all_cell$Ncam1/all_cell$Freq
hist(all_cell$Ncam1_prop, breaks = 50)
#FeaturePlot(data.n, genes_to_ens("Ncam1"), label = T, order = T)
true_Ncam1 <- subset(all_cell, subset = all_cell$Ncam1_prop<=0.52)
write.table(true_Ncam1$Var1, "Myelinated_Ncam1_pc50_res1.4.txt", sep = "\t", quote = F)

true_Ncam1 <- subset(all_cell, subset = all_cell$Ncam1_prop>=0.75)
write.table(true_Ncam1$Var1, "Unmyelinated_Ncam1_high_pc50_res1.4.txt", sep = "\t", quote = F)

##Unmyelinated / lightly myelinated neurons
##Low in Cntnap1, Cntn1; High in Ncam1 // Or low in three of them
hist(all_cell$Cntnap1_prop, breaks = 50)
#FeaturePlot(data.n, genes_to_ens("Cntnap1"), label = T, order = T)
true_Cntnap1 <- subset(all_cell, subset = all_cell$Cntnap1_prop<=0.1)
write.table(true_Cntnap1$Var1, "Unmyelinated_Cntnap1_pc50_res1.4.txt", sep = "\t", quote = F)

hist(all_cell$Cntn1_prop, breaks = 50)
#FeaturePlot(data.n, genes_to_ens("Cntn1"), label = T, order = T)
true_Cntn1 <- subset(all_cell, subset = all_cell$Cntn1_prop<=0.3)
write.table(true_Cntn1$Var1, "Umyelinated_Cntn1_pc50_res1.4.txt", sep = "\t", quote = F)


# -------------------------------------------------------------------------
# Mechanosensor (Piezo2) vs Nocisensor (Trpv1) ---------------------------------------------
##A-beta_LTMRs (Slc17a7, P2ry1, Trpv1(negative))
table(GetAssayData(data.n, slot = "data")[genes_to_ens("Gfra1"),]>0)
Gfra1 <- subset(data.n, subset = ENSMUSG00000025089 >0)
Gfra1_cell <- data.frame(table(Idents(Gfra1)))
all_cell$Gfra1 <- Gfra1_cell$Freq[match(all_cell$Var1, Gfra1_cell$Var1)]
all_cell$Gfra1_prop <- all_cell$Gfra1/all_cell$Freq
hist(all_cell$Gfra1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Gfra1"), label = T, order = T, repel = T)
true_Gfra1 <- subset(all_cell, subset = all_cell$Gfra1_prop>=0.3)
write.table(true_Gfra1$Var1, "Mechano_Gfra1_pc50_res1.4.txt", sep = "\t", quote = F)

table(GetAssayData(data.n, slot = "data")[genes_to_ens("Slc17a7"),]>0)
Slc17a7 <- subset(data.n, subset = ENSMUSG00000070570 >0)
Slc17a7_cell <- data.frame(table(Idents(Slc17a7)))
all_cell$Slc17a7 <- Slc17a7_cell$Freq[match(all_cell$Var1, Slc17a7_cell$Var1)]
all_cell$Slc17a7_prop <- all_cell$Slc17a7/all_cell$Freq
hist(all_cell$Slc17a7_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Slc17a7"), label = T, order = T, repel = T)
true_Slc17a7 <- subset(all_cell, subset = all_cell$Slc17a7_prop>=0.3)
write.table(true_Slc17a7$Var1, "Mechano_Slc17a7_pc50_res1.4.txt", sep = "\t", quote = F)

table(GetAssayData(data.n, slot = "data")[genes_to_ens("Ptgfr"),]>0)
Ptgfr <- subset(data.n, subset = ENSMUSG00000028036 >0)
Ptgfr_cell <- data.frame(table(Idents(Ptgfr)))
all_cell$Ptgfr <- Ptgfr_cell$Freq[match(all_cell$Var1, Ptgfr_cell$Var1)]
all_cell$Ptgfr_prop <- all_cell$Ptgfr/all_cell$Freq
hist(all_cell$Ptgfr_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Ptgfr"), label = T, order = T, repel = T)
true_Ptgfr <- subset(all_cell, subset = all_cell$Ptgfr_prop>=0.1)
write.table(true_Ptgfr$Var1, "Mechano_Ptgfr_pc50_res1.4.txt", sep = "\t", quote = F)

table(GetAssayData(data.n, slot = "data")[genes_to_ens("P2ry1"),]>0)
P2ry1 <- subset(data.n, subset = ENSMUSG00000027765 >0)
P2ry1_cell <- data.frame(table(Idents(P2ry1)))
all_cell$P2ry1 <- P2ry1_cell$Freq[match(all_cell$Var1, P2ry1_cell$Var1)]
all_cell$P2ry1_prop <- all_cell$P2ry1/all_cell$Freq
hist(all_cell$P2ry1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("P2ry1"), label = T, order = T, repel = T)
true_P2ry1 <- subset(all_cell, subset = all_cell$P2ry1_prop>=0.5)
write.table(true_P2ry1$Var1, "Mechano_P2ry1_pc50_res1.4.txt", sep = "\t", quote = F)

table(GetAssayData(data.n, slot = "data")[genes_to_ens("Piezo2"),]>0)
Piezo2 <- subset(data.n, subset = ENSMUSG00000041482 >0)
Piezo2_cell <- data.frame(table(Idents(Piezo2)))
all_cell$Piezo2 <- Piezo2_cell$Freq[match(all_cell$Var1, Piezo2_cell$Var1)]
all_cell$Piezo2_prop <- all_cell$Piezo2/all_cell$Freq
hist(all_cell$Piezo2_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Piezo2"), label = T, order = T, repel = T)
true_Piezo2 <- subset(all_cell, subset = all_cell$Piezo2_prop>=0.5)
write.table(true_Piezo2$Var1, "Mechano_Piezo2_pc50_res1.4.txt", sep = "\t", quote = F)

table(GetAssayData(data.n, slot = "data")[genes_to_ens("Ntrk1"),]>0)
Ntrk1 <- subset(data.n, subset = ENSMUSG00000028072 >0)
Ntrk1_cell <- data.frame(table(Idents(Ntrk1)))
all_cell$Ntrk1 <- Ntrk1_cell$Freq[match(all_cell$Var1, Ntrk1_cell$Var1)]
all_cell$Ntrk1_prop <- all_cell$Ntrk1/all_cell$Freq
hist(all_cell$Ntrk1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Ntrk1"), label = T, order = T, repel = T)
true_Ntrk1 <- subset(all_cell, subset = all_cell$Ntrk1_prop>=0.5)
write.table(true_Ntrk1$Var1, "Nocis_Ntrk1_pc50_res1.4.txt", sep = "\t", quote = F)

table(GetAssayData(data.n, slot = "data")[genes_to_ens("Calca"),]>0)
Calca <- subset(data.n, subset = ENSMUSG00000030669 >0)
Calca_cell <- data.frame(table(Idents(Calca)))
all_cell$Calca <- Calca_cell$Freq[match(all_cell$Var1, Calca_cell$Var1)]
all_cell$Calca_prop <- all_cell$Calca/all_cell$Freq
hist(all_cell$Calca_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Calca"), label = T, order = T, repel = T)
true_Calca <- subset(all_cell, subset = all_cell$Calca_prop>=0.6)
write.table(true_Calca$Var1, "Nocis_Calca_pc50_res1.4.txt", sep = "\t", quote = F)

table(GetAssayData(data.n, slot = "data")[genes_to_ens("Trpv1"),]>0)
Trpv1 <- subset(data.n, subset = ENSMUSG00000005952 >0)
Trpv1_cell <- data.frame(table(Idents(Trpv1)))
all_cell$Trpv1 <- Trpv1_cell$Freq[match(all_cell$Var1, Trpv1_cell$Var1)]
all_cell$Trpv1_prop <- all_cell$Trpv1/all_cell$Freq
hist(all_cell$Trpv1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Trpv1"), label = T, order = T, repel = T)
true_Trpv1 <- subset(all_cell, subset = all_cell$Trpv1_prop>=0.4)
write.table(true_Trpv1$Var1, "Nocis_Trpv1_pc50_res1.4.txt", sep = "\t", quote = F)


# Organ innervated (Zhao) -------------------------------------------------
###Select the clusters with co-expression
Sox4 <- subset(data.n, subset = ENSMUSG00000076431 >0)
Sox4_cell <- data.frame(table(Idents(Sox4)))
all_cell$Sox4 <- Sox4_cell$Freq[match(all_cell$Var1, Sox4_cell$Var1)]
all_cell$Sox4_prop <- all_cell$Sox4/all_cell$Freq
hist(all_cell$Sox4_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Sox4"), label = T, order = T)
true_Sox4 <- subset(all_cell, subset = all_cell$Sox4_prop>=0.5)
write.table(true_Sox4$Var1, "Lung_Sox4_pc50_res1.4.txt", sep = "\t", quote = F)

Runx1 <- subset(data.n, subset = ENSMUSG00000022952 >0)
Runx1_cell <- data.frame(table(Idents(Runx1)))
all_cell$Runx1 <- Runx1_cell$Freq[match(all_cell$Var1, Runx1_cell$Var1)]
all_cell$Runx1_prop <- all_cell$Runx1/all_cell$Freq
hist(all_cell$Runx1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Runx1"), label = T, order = T)
true_Runx1 <- subset(all_cell, subset = all_cell$Runx1_prop>=0.6)
write.table(true_Runx1$Var1, "Lung_Runx1_pc50_res1.4.txt", sep = "\t", quote = F)

Pou4f1 <- subset(data.n, subset = ENSMUSG00000048349 >0)
Pou4f1_cell <- data.frame(table(Idents(Pou4f1)))
all_cell$Pou4f1 <- Pou4f1_cell$Freq[match(all_cell$Var1, Pou4f1_cell$Var1)]
all_cell$Pou4f1_prop <- all_cell$Pou4f1/all_cell$Freq
hist(all_cell$Pou4f1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Pou4f1"), label = T, order = T)
true_Pou4f1 <- subset(all_cell, subset = all_cell$Pou4f1_prop>=0.7)
write.table(true_Pou4f1$Var1, "Lung_Pou4f1_pc50_res1.4.txt", sep = "\t", quote = F)

####
Irf6 <- subset(data.n, subset = ENSMUSG00000026638 >0)
Irf6_cell <- data.frame(table(Idents(Irf6)))
all_cell$Irf6 <- Irf6_cell$Freq[match(all_cell$Var1, Irf6_cell$Var1)]
all_cell$Irf6_prop <- all_cell$Irf6/all_cell$Freq
hist(all_cell$Irf6_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Irf6"), label = T, order = T)
true_Irf6 <- subset(all_cell, subset = all_cell$Irf6_prop>=0.2)
write.table(true_Irf6$Var1, "Heart_Irf6_pc50_res1.4.txt", sep = "\t", quote = F)

Esr1 <- subset(data.n, subset = ENSMUSG00000019768 >0)
Esr1_cell <- data.frame(table(Idents(Esr1)))
all_cell$Esr1 <- Esr1_cell$Freq[match(all_cell$Var1, Esr1_cell$Var1)]
all_cell$Esr1_prop <- all_cell$Esr1/all_cell$Freq
hist(all_cell$Esr1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Esr1"), label = T, order = T)
true_Esr1 <- subset(all_cell, subset = all_cell$Esr1_prop>=0.35)
write.table(true_Esr1$Var1, "Heart_Esr1_pc50_res1.4.txt", sep = "\t", quote = F)

Tbx3 <- subset(data.n, subset = ENSMUSG00000018604 >0)
Tbx3_cell <- data.frame(table(Idents(Tbx3)))
all_cell$Tbx3 <- Tbx3_cell$Freq[match(all_cell$Var1, Tbx3_cell$Var1)]
all_cell$Tbx3_prop <- all_cell$Tbx3/all_cell$Freq
hist(all_cell$Tbx3_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Tbx3"), label = T, order = T)
true_Tbx3 <- subset(all_cell, subset = all_cell$Tbx3_prop>=0.6)
write.table(true_Tbx3$Var1, "Heart_Tbx3_pc50_res1.4.txt", sep = "\t", quote = F)

Mef2c <- subset(data.n, subset = ENSMUSG00000005583 >0)
Mef2c_cell <- data.frame(table(Idents(Mef2c)))
all_cell$Mef2c <- Mef2c_cell$Freq[match(all_cell$Var1, Mef2c_cell$Var1)]
all_cell$Mef2c_prop <- all_cell$Mef2c/all_cell$Freq
hist(all_cell$Mef2c_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Mef2c"), label = T, order = T)
true_Mef2c <- subset(all_cell, subset = all_cell$Mef2c_prop>=0.4)
write.table(true_Mef2c$Var1, "Heart_Mef2c_pc50_res1.4.txt", sep = "\t", quote = F)

###
Etv1 <- subset(data.n, subset = ENSMUSG00000004151 >0)
Etv1_cell <- data.frame(table(Idents(Etv1)))
all_cell$Etv1 <- Etv1_cell$Freq[match(all_cell$Var1, Etv1_cell$Var1)]
all_cell$Etv1_prop <- all_cell$Etv1/all_cell$Freq
hist(all_cell$Etv1_prop, breaks = 50)
#FeaturePlot(data.n, genes_to_ens("Etv1"), label = T, order = T)
true_Etv1 <- subset(all_cell, subset = all_cell$Etv1_prop>=0.8)
write.table(true_Etv1$Var1, "Gut_Etv1_pc50_res1.4.txt", sep = "\t", quote = F)

Zfhx3 <- subset(data.n, subset = ENSMUSG00000038872 >0)
Zfhx3_cell <- data.frame(table(Idents(Zfhx3)))
all_cell$Zfhx3 <- Zfhx3_cell$Freq[match(all_cell$Var1, Zfhx3_cell$Var1)]
all_cell$Zfhx3_prop <- all_cell$Zfhx3/all_cell$Freq
hist(all_cell$Zfhx3_prop, breaks = 50)
#FeaturePlot(data.n, genes_to_ens("Zfhx3"), label = T, order = T)
true_Zfhx3 <- subset(all_cell, subset = all_cell$Zfhx3_prop>=0.8)
write.table(true_Zfhx3$Var1, "Gut_Zfhx3_pc50_res1.4.txt", sep = "\t", quote = F)

Rab15 <- subset(data.n, subset = ENSMUSG00000021062 >0)
Rab15_cell <- data.frame(table(Idents(Rab15)))
all_cell$Rab15 <- Rab15_cell$Freq[match(all_cell$Var1, Rab15_cell$Var1)]
all_cell$Rab15_prop <- all_cell$Rab15/all_cell$Freq
hist(all_cell$Rab15_prop, breaks = 50)
#FeaturePlot(data.n, genes_to_ens("Rab15"), label = T, order = T)
true_Rab15 <- subset(all_cell, subset = all_cell$Rab15_prop>=0.7)
write.table(true_Rab15$Var1, "Gut_Rab15_pc50_res1.4.txt", sep = "\t", quote = F)

###
Nhlh2 <- subset(data.n, subset = ENSMUSG00000048540 >0)
Nhlh2_cell <- data.frame(table(Idents(Nhlh2)))
all_cell$Nhlh2 <- Nhlh2_cell$Freq[match(all_cell$Var1, Nhlh2_cell$Var1)]
all_cell$Nhlh2_prop <- all_cell$Nhlh2/all_cell$Freq
hist(all_cell$Nhlh2_prop, breaks = 50)
#FeaturePlot(data.n, genes_to_ens("Nhlh2"), label = T, order = T)
true_Nhlh2 <- subset(all_cell, subset = all_cell$Nhlh2_prop>=0.4)
write.table(true_Nhlh2$Var1, "Pancreas_Nhlh2_pc50_res1.4.txt", sep = "\t", quote = F)

Klf4 <- subset(data.n, subset = ENSMUSG00000003032 >0)
Klf4_cell <- data.frame(table(Idents(Klf4)))
all_cell$Klf4 <- Klf4_cell$Freq[match(all_cell$Var1, Klf4_cell$Var1)]
all_cell$Klf4_prop <- all_cell$Klf4/all_cell$Freq
hist(all_cell$Klf4_prop, breaks = 50)
#FeaturePlot(data.n, genes_to_ens("Klf4"), label = T, order = T)
true_Klf4 <- subset(all_cell, subset = all_cell$Klf4_prop>=0.15)
write.table(true_Klf4$Var1, "Pancreas_Klf4_pc50_res1.4.txt", sep = "\t", quote = F)

#Gut innervated (Bai) ------------------------------------------------
table(GetAssayData(data.n, slot = "data")[genes_to_ens("Sst"),]>0)
Sst <- subset(data.n, subset = ENSMUSG00000004366 >0)
Sst_cell <- data.frame(table(Idents(Sst)))
all_cell$Sst <- Sst_cell$Freq[match(all_cell$Var1, Sst_cell$Var1)]
all_cell$Sst_prop <- all_cell$Sst/all_cell$Freq
hist(all_cell$Sst_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Sst"), label = T, order = T)
true_Sst <- subset(all_cell, subset = all_cell$Sst_prop>=0.1)
write.table(true_Sst$Var1, "Stomach_Sst_pc50_res1.4.txt", sep = "\t", quote = F)

table(GetAssayData(data.n, slot = "data")[genes_to_ens("Calca"),]>0)
Calca <- subset(data.n, subset = ENSMUSG00000030669 >0)
Calca_cell <- data.frame(table(Idents(Calca)))
all_cell$Calca <- Calca_cell$Freq[match(all_cell$Var1, Calca_cell$Var1)]
all_cell$Calca_prop <- all_cell$Calca/all_cell$Freq
hist(all_cell$Calca_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Calca"), label = T, order = T)
true_Calca <- subset(all_cell, subset = all_cell$Calca_prop>=0.6)
write.table(true_Calca$Var1, "Stomach_Calca_pc50_res1.4.txt", sep = "\t", quote = F)

###
Dbh <- subset(data.n, subset = ENSMUSG00000000889 >0)
Dbh_cell <- data.frame(table(Idents(Dbh)))
all_cell$Dbh <- Dbh_cell$Freq[match(all_cell$Var1, Dbh_cell$Var1)]
all_cell$Dbh_prop <- all_cell$Dbh/all_cell$Freq
hist(all_cell$Dbh_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Dbh"), label = T, order = T)
true_Dbh <- subset(all_cell, subset = all_cell$Dbh_prop>=0.2)
write.table(true_Dbh$Var1, "Duodenum_Dbh_pc50_res1.4.txt", sep = "\t", quote = F)

Ebf3 <- subset(data.n, subset = ENSMUSG00000010476 >0)
Ebf3_cell <- data.frame(table(Idents(Ebf3)))
all_cell$Ebf3 <- Ebf3_cell$Freq[match(all_cell$Var1, Ebf3_cell$Var1)]
all_cell$Ebf3_prop <- all_cell$Ebf3/all_cell$Freq
hist(all_cell$Ebf3_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Ebf3"), label = T, order = T)
true_Ebf3 <- subset(all_cell, subset = all_cell$Ebf3_prop>=0.2)
write.table(true_Ebf3$Var1, "Duodenum_Ebf3_pc50_res1.4.txt", sep = "\t", quote = F)

###
Caln1 <- subset(data.n, subset = ENSMUSG00000060371 >0)
Caln1_cell <- data.frame(table(Idents(Caln1)))
all_cell$Caln1 <- Caln1_cell$Freq[match(all_cell$Var1, Caln1_cell$Var1)]
all_cell$Caln1_prop <- all_cell$Caln1/all_cell$Freq
hist(all_cell$Caln1_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Caln1"), label = T, order = T)
true_Caln1 <- subset(all_cell, subset = all_cell$Caln1_prop>=0.05)
write.table(true_Caln1$Var1, "Jejunum_Ileum_Caln1_pc50_res1.4.txt", sep = "\t", quote = F)

Edn3 <- subset(data.n, subset = ENSMUSG00000027524 >0)
Edn3_cell <- data.frame(table(Idents(Edn3)))
all_cell$Edn3 <- Edn3_cell$Freq[match(all_cell$Var1, Edn3_cell$Var1)]
all_cell$Edn3_prop <- all_cell$Edn3/all_cell$Freq
hist(all_cell$Edn3_prop, breaks = 20)
#FeaturePlot(data.n, genes_to_ens("Edn3"), label = T, order = T)
true_Edn3 <- subset(all_cell, subset = all_cell$Edn3_prop>=0.05)
write.table(true_Edn3$Var1, "Jejunum_Ileum_Edn3_pc50_res1.4.txt", sep = "\t", quote = F)


###
Areg <- subset(data.n, subset = ENSMUSG00000029378 >0)
Areg_cell <- data.frame(table(Idents(Areg)))
all_cell$Areg <- Areg_cell$Freq[match(all_cell$Var1, Areg_cell$Var1)]
all_cell$Areg_prop <- all_cell$Areg/all_cell$Freq
hist(all_cell$Areg_prop, breaks = 20)
true_Areg <- subset(all_cell, subset = all_cell$Areg_prop>=0.1)
write.table(true_Areg$Var1, "Colon_Areg_pc50_res1.4.txt", sep = "\t", quote = F)



###Sodium channel
sodium_channel <- data.frame(row.names = all_cell$Var1,
                             Scn1a = all_cell$Scn1a_prop,
                             Scn10a = all_cell$Scn10a_prop)

#Add classification
classification <- data.frame(row.names = all_cell$Var1,
                             Sodium_channel = c("Nav1.8",
                                                "Nav1.8",
                                                "Nav1.1/Nav1.8",
                                                "Nav1.8",
                                                "Nav1.1",
                                                "Nav1.8",
                                                "Nav1.8",
                                                "Nav1.8",
                                                "Nav1.8",
                                                "Nav1.8",
                                                "Nav1.1",
                                                "Nav1.1",
                                                "Nav1.8",
                                                "Nav1.8",
                                                "Nav1.8",
                                                "Nav1.8",
                                                "Nav1.1",
                                                "Nav1.8",
                                                "Nav1.1",
                                                "Nav1.1",
                                                "Nav1.1",
                                                "Nav1.1",
                                                "Nav1.1/Nav1.8",
                                                "Nav1.1",
                                                "Nav1.8",
                                                "Nav1.1"))

#Set classification colour
col_class <- list(Sodium_channel = brewer.pal(8, "Set2")[2:4])
names(col_class$Sodium_channel) <- unique(classification$Sodium_channel)

pdf("nav_heatmap.pdf", width = 10, height = 15)
pheatmap(sodium_channel,color = colorRampPalette(rev(brewer.pal(n=9, name = "PRGn")))(100),
         cluster_rows = T,
         cluster_cols = F,
         angle_col = 315,
         cellwidth = 25,
         cellheight = 25,
         annotation_row = classification,
         annotation_colors = col_class,
         fontfamily = "Helvetica",
         fontsize = 20)
dev.off()

pdf("nav_dimplot.pdf", width = 10, height = 10)
DimPlot(data.n, group.by = "Sodium_channel", raster = F, order = T,
        cols = col_class$Sodium_channel) + 
  ggtitle('') + NoAxes() + NoLegend() + coord_fixed()
dev.off()




###Fibre type
fibre_type <- data.frame(row.names = all_cell$Var1,
                         Cntnap1 = all_cell$Cntnap1_prop,
                         Cntn1 = all_cell$Cntn1_prop,
                         Nefh = all_cell$Nefh_prop,
                         Ncam1 = all_cell$Ncam1_prop)

classification <- data.frame(row.names = all_cell$Var1,
                             Fibre_type = c("Unmyelinated",
                                            "Unmyelinated",
                                            "Myelinated",
                                            "Lightly myelinated",
                                            "Lightly myelinated",
                                            "Lightly myelinated",
                                            "Lightly myelinated",
                                            "Unmyelinated",
                                            "Unmyelinated",
                                            "Unmyelinated",
                                            "Lightly myelinated",
                                            "Lightly myelinated",
                                            "Unmyelinated",
                                            "Unmyelinated",
                                            "Myelinated",
                                            "Unmyelinated",
                                            "Myelinated",
                                            "Unmyelinated",
                                            "Myelinated",
                                            "Lightly myelinated",
                                            "Lightly myelinated",
                                            "Lightly myelinated",
                                            "Lightly myelinated",
                                            "Unmyelinated",
                                            "Unmyelinated",
                                            "Myelinated"))

col_class <- list(Fibre_type = brewer.pal(8, "Set1")[1:3])
names(col_class$Fibre_type) <- unique(classification$Fibre_type)

pdf("fibre_heatmap.pdf", width = 10, height = 15)
pheatmap(fibre_type, color = colorRampPalette(rev(brewer.pal(n=9, name = "PuOr")))(100),
         cluster_rows = T,
         cluster_cols = F,
         angle_col = 315,
         cellwidth = 25,
         cellheight = 25,
         annotation_row = classification,
         annotation_colors = col_class,
         fontfamily = "Helvetica",
         fontsize = 20)
dev.off()

pdf("fibre_dimplot.pdf", width = 10, height = 10)
DimPlot(data.n, group.by = "Fibre_type", raster = F, order = T,
        cols = col_class$Fibre_type) + 
  ggtitle('') + NoAxes() + NoLegend() + coord_fixed()
dev.off()

###Sensor type
sensor_type <- data.frame(row.names = all_cell$Var1,
                          Gfra1 = all_cell$Gfra1_prop,
                          Slc17a7 = all_cell$Slc17a7_prop,
                          Ptgfr = all_cell$Ptgfr_prop,
                          P2ry1 = all_cell$P2ry1_prop, 
                          Piezo2 = all_cell$Piezo2_prop,
                          Ntrk1 = all_cell$Ntrk1_prop,
                          Calca = all_cell$Calca_prop, 
                          Trpv1 = all_cell$Trpv1_prop)

classification <- data.frame(row.names = all_cell$Var1,
                             Sensor_type = c("Mechanosensor/Nocisensor",
                                             "Nocisensor",
                                             "Mechanosensor/Nocisensor",
                                             "Mechanosensor",
                                             "Nocisensor",
                                             "Nocisensor",
                                             "Nocisensor",
                                             "Mechanosensor",
                                             "Nocisensor",
                                             "Mechanosensor/Nocisensor",
                                             "Mechanosensor",
                                             "Mechanosensor",
                                             "Mechanosensor/Nocisensor",
                                             "Nocisensor",
                                             "Mechanosensor",
                                             "Nocisensor",
                                             "Mechanosensor/Nocisensor",
                                             "Nocisensor",
                                             "Mechanosensor",
                                             "Mechanosensor",
                                             "Mechanosensor",
                                             "Mechanosensor",
                                             "Mechanosensor/Nocisensor",
                                             "Nocisensor",
                                             "Nocisensor",
                                             "Mechanosensor"))

col_class <- list(Sensor_type = brewer.pal(8, "Dark2")[1:3])
names(col_class$Sensor_type) <- unique(classification$Sensor_type)

pdf("sensor_heatmap.pdf", width = 10, height = 15)
pheatmap(sensor_type, color = colorRampPalette(rev(brewer.pal(n=9, name = "Spectral")))(100),
         cluster_rows = T,
         cluster_cols = F,
         angle_col = 315,
         cellwidth = 25,
         cellheight = 25,
         annotation_row = classification,
         annotation_colors = col_class,
         fontfamily = "Helvetica",
         fontsize = 20)
dev.off()

pdf("sensor_dimplot.pdf", width = 10, height = 10)
DimPlot(data.n, group.by = "Sensor_type", raster = F, order = T,
        cols = col_class$Sensor_type) + 
  ggtitle('') + NoAxes() + NoLegend() + coord_fixed()
dev.off()


###Organ projection
organ_projection <- data.frame(row.names = all_cell$Var1,
                               Sox4 = all_cell$Sox4_prop,
                               Runx1 = all_cell$Runx1_prop,
                               Pou4f1 = all_cell$Pou4f1_prop,
                               Irf6 = all_cell$Irf6_prop,
                               Esr1 = all_cell$Esr1_prop,
                               Tbx3 = all_cell$Tbx3_prop,
                               Mef2c = all_cell$Mef2c_prop,
                               Nhlh2 = all_cell$Nhlh2_prop,
                               Klf4 = all_cell$Klf4_prop,
                               Etv1 = all_cell$Etv1_prop,
                               Zfhx3 = all_cell$Zfhx3_prop,
                               Rab15 = all_cell$Rab15_prop,
                               Sst = all_cell$Sst_prop,
                               Calca = all_cell$Calca_prop,
                               Dbh = all_cell$Dbh_prop,
                               Ebf3 = all_cell$Ebf3_prop,
                               Caln1 = all_cell$Caln1_prop,
                               Edn3 = all_cell$Edn3_prop,
                               Areg = all_cell$Areg_prop)

classification <- data.frame(row.names = all_cell$Var1,
                             Organ_projection = c("Broad projection",
                                                  "Broad projection",
                                                  "Gut",
                                                  "Broad projection",
                                                  "Lung",
                                                  "Broad projection",
                                                  "Broad projection",
                                                  "Gut",
                                                  "Broad projection",
                                                  "Broad projection",
                                                  "Broad projection",
                                                  "Pancreas",
                                                  "Broad projection",
                                                  "Gut",
                                                  "Gut",
                                                  "Gut",
                                                  "Duodenum",
                                                  "Broad projection",
                                                  "Gut",
                                                  "Heart",
                                                  "Gut",
                                                  "Broad projection",
                                                  "Gut",
                                                  "Gut",
                                                  "Broad projection",
                                                  "Jejunum/Ileum"))

col_class <- list(Organ_projection = brewer.pal(7, "Accent")[1:7])
names(col_class$Organ_projection) <- unique(classification$Organ_projection)

pdf("organ_heatmap.pdf", width = 15, height = 15)
pheatmap(organ_projection, color = colorRampPalette(rev(brewer.pal(n=9, name = "RdYlBu")))(100),
         cluster_rows = T,
         cluster_cols = F,
         angle_col = 315,
         cellwidth = 25,
         cellheight = 25,
         annotation_row = classification,
         annotation_colors = col_class,
         fontfamily = "Helvetica",
         fontsize = 20)
dev.off()

pdf("organ_dimplot.pdf", width = 10, height = 10)
DimPlot(data.n, group.by = "Organ_projection", raster = F, order = T,
        cols = col_class$Organ_projection) + 
  ggtitle('') + NoAxes() + NoLegend() + coord_fixed()
dev.off()




