library(Seurat)
library(dplyr)

#Use the markers identified by Wilcox (Seurat_Findallmarkers)
markers <- read.delim("251019_nodose_markers.txt")
markers$prop <- markers$avg_log2FC*(markers$pct.1 + 0.01)/(markers$pct.2 + 0.01)
markers$genesymbol <- gsub("\\..*","", markers$genesymbol)

top3 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = prop)

write.table(markers, file = "Wilcox_updated_2510.txt", sep = "\t", quote = F, row.names = F)

clustermarkers <- c("EC1.Cldn5/Ly6c1/Ly6a",
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
