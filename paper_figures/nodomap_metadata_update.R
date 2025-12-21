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


# change data idents to cluster_number and order ----------------------------------------------

data <- SetIdent(data, value = "cluster_number")
Idents(data) <- factor(Idents(data), levels = celltypes)
levels(Idents(data))
sodium_channel <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                    NA,NA,NA,NA,NA,NA,
                    "Nav1.8",
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
                    "Nav1.1")
names(sodium_channel) <- levels(data)
data <- RenameIdents(data,sodium_channel)
data <- AddMetaData(data, metadata = Idents(data), col.name = "Sodium_channel")


data <- SetIdent(data, value = "cluster_number")
Idents(data) <- factor(Idents(data), levels = celltypes)
levels(Idents(data))
fibre_type <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                NA,NA,NA,NA,NA,NA,
                "Unmyelinated",
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
                "Myelinated")
names(fibre_type) <- levels(data)
data <- RenameIdents(data,fibre_type)
data <- AddMetaData(data, metadata = Idents(data), col.name = "Fibre_type")


data <- SetIdent(data, value = "cluster_number")
Idents(data) <- factor(Idents(data), levels = celltypes)
levels(Idents(data))
sensor_type <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                 NA,NA,NA,NA,NA,NA,
                 "Mechanosensor/Nocisensor",
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
                 "Mechanosensor")
names(sensor_type) <- levels(data)
data <- RenameIdents(data,sensor_type)
data <- AddMetaData(data, metadata = Idents(data), col.name = "Sensor_type")


data <- SetIdent(data, value = "cluster_number")
Idents(data) <- factor(Idents(data), levels = celltypes)
levels(Idents(data))
organ_innervated <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                      NA,NA,NA,NA,NA,NA,
                      "Broad projection",
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
                      "Jejunum/Ileum")
names(organ_innervated) <- levels(data)
data <- RenameIdents(data,organ_innervated)
data <- AddMetaData(data, metadata = Idents(data), col.name = "Organ_projection")

View(data@meta.data)
saveRDS(data, "diet_nodomap_metadata_updated_2510.rds")


#RCTD
data.st <- readRDS("nodose_st_RCTD_paramrelax_251025.RDS")
DefaultAssay(data.st) <- "RNA"
data.st <- SetIdent(data.st, value = "first_type") 
Idents(data.st) <- factor(Idents(data.st), levels = celltypes)
levels(Idents(data.st))

x<-data.st@meta.data[,c(17,18)]
levels(x$first_type)
x$first_type[x$spot_class != 'singlet'] <- NA
colnames(x)<-c('spot_class', 'singlet_cells')
data.st <- AddMetaData(data.st, metadata = x$singlet_cells, col.name = "singlet_cells")

data.st <- SetIdent(data.st, value = "singlet_cells") 
Idents(data.st) <- factor(Idents(data.st), levels = celltypes)
levels(Idents(data.st))

sodium_channel <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                    NA,NA,NA,NA,NA,NA,NA,NA,NA,
                    "Nav1.8",
                    "Nav1.8",
                    "Nav1.1/Nav1.8",
                    "Nav1.8",
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
                    "Nav1.1")
names(sodium_channel) <- levels(data.st)
data.st <- RenameIdents(data.st,sodium_channel)
data.st <- AddMetaData(data.st, metadata = Idents(data.st), col.name = "Sodium_channel")


data.st <- SetIdent(data.st, value = "singlet_cells") 
Idents(data.st) <- factor(Idents(data.st), levels = celltypes)
levels(Idents(data.st))
fibre_type <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                NA,NA,NA,NA,NA,NA,NA,NA,NA,
                "Unmyelinated",
                "Unmyelinated",
                "Myelinated",
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
                "Myelinated")
names(fibre_type) <- levels(data.st)
data.st <- RenameIdents(data.st,fibre_type)
data.st <- AddMetaData(data.st, metadata = Idents(data.st), col.name = "Fibre_type")


data.st <- SetIdent(data.st, value = "singlet_cells") 
Idents(data.st) <- factor(Idents(data.st), levels = celltypes)
levels(Idents(data.st))
sensor_type <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                 NA,NA,NA,NA,NA,NA,NA,NA,NA,
                 "Mechanosensor/Nocisensor",
                 "Nocisensor",
                 "Mechanosensor/Nocisensor",
                 "Mechanosensor",
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
                 "Mechanosensor")
names(sensor_type) <- levels(data.st)
data.st <- RenameIdents(data.st,sensor_type)
data.st <- AddMetaData(data.st, metadata = Idents(data.st), col.name = "Sensor_type")


data.st <- SetIdent(data.st, value = "singlet_cells") 
Idents(data.st) <- factor(Idents(data.st), levels = celltypes)
levels(Idents(data.st))
organ_innervated <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
                      NA,NA,NA,NA,NA,NA,NA,NA,NA,
                      "Broad projection",
                      "Broad projection",
                      "Gut",
                      "Broad projection",
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
                      "Jejunum/Ileum")

names(organ_innervated) <- levels(data.st)
data.st <- RenameIdents(data.st,organ_innervated)
data.st <- AddMetaData(data.st, metadata = Idents(data.st), col.name = "Organ_projection")

View(data.st@meta.data)
saveRDS(data.st, "nodose_st_RCTD_251025_metadata_updated.RDS")











