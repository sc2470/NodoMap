library(Seurat)
library(dplyr)

#Use the markers identified by Wilcox (Seurat_Findallmarkers)
markers <- read.delim("pc50-resolution-1.4-markers-230613.txt")
markers$prop <- markers$avg_log2FC*(markers$pct.1 + 0.01)/(markers$pct.2 + 0.01)
markers$genesymbol <- gsub("\\..*","", markers$genesymbol)

top3 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = prop)

write.table(markers, file = "pc50-resolution-1.4-markers-Wilcox.txt", sep = "\t", quote = F, row.names = F)

#Load the data
data <- readRDS("../nodose_integrated_PC30-60_res0.6-2/nodose-integrated-pcs-50-resolution-1.4.RDS")
DefaultAssay(data) <- "originalexp"
data <- SetIdent(data, value = "integrated_snn_res.1.4")

mouse_gene_names<-read.delim("../../mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

#Cluster markers
clustermarkers <- c("EC1.Cldn5/Ly6c1/Ly6a",
                    "SGC1.Mmd2/Ednrb/Fabp7",
                    "SGC2.Cxcl10/Acsbg1/Phgdh",
                    "MGC1.Ncmap/Fam178b/Cldn19",
                    "SGC3.Bcan/Fbln2/Fbln5",
                    "NGN1.C1ql3/Adcyap1/St18",
                    "HC1.C1qc/C1qa/C1qb",
                    "SGC4.Atf3/Emp1/Cebpd",
                    "SGC5.Pou3f1/Col1a2/Gas7",
                    "FB1.Lum/Ccl11/Smoc2",
                    "SGC6.Ifit3/Cdh19/Isg15",
                    "MGC2.Prx/Pllp/Bcas1",
                    "FB2.Acta2/Myh11/Rgs5",
                    "NGN2.Htr3b/Htr3a/Hs3st4",
                    "GC1.Kcna1/Fxyd3/Mbp",
                    "NGN3.Cysltr2/Chrnb3/Kcnip4",
                    "MGC3.Slc36a2/Nr4a2/Ugt8a",
                    "NGN4.Uts2b/Cckar/Vip",
                    "SGC7.G0s2/Sdc4/Ccn1",
                    "MGC4.Mt2/Tcim/Drp2",
                    "NGN5.Rtn1/Ngfr/Clu",
                    "NGN6.Kcng1/Vmn1r85/Trpv1",
                    "SGC8.Ttyh1/Ptprz1/Sfrp5",
                    "JGN1.Mrgprd/Tmem45b/Grik1",
                    "NGN7.Lypd1/Ddc/Tafa2",
                    "GC2.Sostdc1/Plekhb1/Tsc22d4",
                    "JGN2.Gfra3/Tac1/Tmem255a",
                    "FB3.Igfbp6/Thbs4/Islr",
                    "NGN8.Rbp4/Gda/Kcnk9",
                    "NGN9.Trpa1/Nos1/Efcab6",
                    "EC2.Selp/Fabp4/Vwf",
                    "NGN10.Miat/Snhg11/Meg3",
                    "NGN11.Olfm3/Tmem233/P2ry1",
                    "NGN12.Sprr1a/Ecel1/Nts",
                    "NGN13.Slc18a3/Pappa2/Runx3",
                    "EC3.Podxl/Depp1/Ptprb",
                    "NGN14.Lox/Pkib/Kcnk3",
                    "NGN15.Slc17a7/Hapln4/Lgi3",
                    "JGN3.Trappc3l/Nptx1/Tafa1",
                    "HC2.Cd52/Lsp1/Rac2",
                    "NGN16.Bmp3/Gal/Sost",
                    "JGN4.Wfdc2/C1ql4/Rarres1",
                    "NGN17.Gabra1/Chodl/Gabrb2",
                    "NGN18.Thsd7b/Cacng5/Gata3",
                    "FB4.Hbb-bs/Bnc2/Ebf2",
                    "GC3.C1qc/Tyrobp/Lyz2",
                    "NGN19.Gabra1/Gabrb2/Chrm2",
                    "JGN5.Trpm8/Foxp2/Gsg1l",
                    "MGC5.Smoc2/Ccl11/Lum",
                    "MGC6.Ifitm1/Kcnj8/Rgs5",
                    "NGN20.Slc6a2/Npy/Hand2",
                    "NGN21.Amigo2/Glp1r/Npy2r",
                    "NGN22.Olfr78/Runx3/Avpr1a")







