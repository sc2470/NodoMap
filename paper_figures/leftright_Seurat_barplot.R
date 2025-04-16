library(ggplot2)
library(reshape2)
library(dplyr)

seurat <- list.files(path = ".", pattern = "^leftright\\d+.txt")
seurat <- do.call(rbind, Map("cbind", lapply(seurat, read.delim), cluster=seurat))
seurat$cluster <- gsub("\\D","",seurat$cluster)
colnames(seurat)[1] <- "gene"

seurat <- subset(seurat, subset = p_val < 0.05)

degs <- seurat %>%
  group_by(cluster) %>%
  dplyr::count(cluster)
up <- seurat %>%
  group_by(cluster) %>%
  summarise(upregulated = sum(avg_log2FC>0))
down <- seurat %>%
  group_by(cluster) %>%
  summarise(downregulated = sum(avg_log2FC<0))
degs$Right_Up <- up$upregulated[match(degs$cluster, up$cluster)]
degs$Right_Down <- down$downregulated[match(degs$cluster, down$cluster)]

degs$cluster <- as.numeric(degs$cluster)
degs <- degs %>%
  arrange(cluster)

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

degs$clustermarkers <- clustermarkers
degs$Right_Down<- degs$Right_Down*-1

temp <- melt(degs, "clustermarkers")
temp <- subset(temp, variable %in% c("Right_Up","Right_Down"))
temp$clustermarkers <- factor(temp$clustermarkers, levels = unique(temp$clustermarkers))

split.row <- gsub("\\d","",temp$clustermarkers)
temp <- split.data.frame(temp, f = factor(split.row))

temp1 <- within(temp, rm("NGN", "JGN"))
temp1 <- do.call(rbind,temp1)
temp1$clustermarkers <- factor(temp1$clustermarkers, levels = unique(temp1$clustermarkers))

pdf("figure_Seurat_DEGs.pdf", width = 7, height = 9)
ggplot(temp$NGN, aes(x = value, y = clustermarkers, fill = variable)) +
  geom_col(stat="stack") + theme_classic() + 
  scale_fill_manual(values = c("Right_Down" = "#4477AA", "Right_Up" = "#AA3377"), 
                    labels = c("Enriched in Left", "Enriched in Right")) +
  xlab("Number of DEGs") +
  ylab(NULL) + guides(fill= guide_legend(byrow = T, title = NULL)) +
  theme(axis.text = element_text(size = 19, color = "black", family = "Helvetica"),
        axis.title = element_text(size = 19),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 19),
        legend.spacing.y = unit(10, 'pt'))

ggplot(temp$JGN, aes(x = value, y = clustermarkers, fill = variable)) +
  geom_col(stat="stack", width = 0.2) + theme_classic() + 
  scale_fill_manual(values = c("Right_Down" = "#4477AA", "Right_Up" = "#AA3377"),
                    labels = c("Enriched in Left", "Enriched in Right")) +
  xlab("Number of DEGs") +
  ylab(NULL) + guides(fill= guide_legend(byrow = T, title = NULL)) +
  theme(axis.text = element_text(size = 19, color = "black", family = "Helvetica"),
        axis.title = element_text(size = 19),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 19),
        legend.spacing.y = unit(10, 'pt'))

ggplot(temp1, aes(x = value, y = clustermarkers, fill = variable)) +
  geom_col(stat="stack") + theme_classic() + 
  scale_fill_manual(values = c("Right_Down" = "#4477AA", "Right_Up" = "#AA3377"),
                    labels = c("Enriched in Left", "Enriched in Right")) +
  xlab("Number of DEGs") +
  ylab(NULL) + guides(fill= guide_legend(byrow = T, title = NULL)) +
  theme(axis.text = element_text(size = 19, color = "black", family = "Helvetica"),
        axis.title = element_text(size = 19),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 19),
        legend.spacing.y = unit(10, 'pt'))
dev.off()

####
seurat$clustermarkers <- clustermarkers[as.integer(seurat$cluster) + 1]
genemarkers <- c("Cldn5/Ly6c1/Ly6a",
                 "Mmd2/Ednrb/Fabp7",
                 "Cxcl10/Acsbg1/Phgdh",
                 "Ncmap/Fam178b/Cldn19",
                 "Bcan/Fbln2/Fbln5",
                 "C1ql3/Adcyap1/St18",
                 "C1qc/C1qa/C1qb",
                 "Atf3/Emp1/Cebpd",
                 "Pou3f1/Col1a2/Gas7",
                 "Lum/Ccl11/Smoc2",
                 "Ifit3/Cdh19/Isg15",
                 "Prx/Pllp/Bcas1",
                 "Acta2/Myh11/Rgs5",
                 "Htr3b/Htr3a/Hs3st4",
                 "Kcna1/Fxyd3/Mbp",
                 "Cysltr2/Chrnb3/Kcnip4",
                 "Slc36a2/Nr4a2/Ugt8a",
                 "Uts2b/Cckar/Vip",
                 "G0s2/Sdc4/Ccn1",
                 "Mt2/Tcim/Drp2",
                 "Rtn1/Ngfr/Clu",
                 "Kcng1/Vmn1r85/Trpv1",
                 "Ttyh1/Ptprz1/Sfrp5",
                 "Mrgprd/Tmem45b/Grik1",
                 "Lypd1/Ddc/Tafa2",
                 "Sostdc1/Plekhb1/Tsc22d4",
                 "Gfra3/Tac1/Tmem255a",
                 "Igfbp6/Thbs4/Islr",
                 "Rbp4/Gda/Kcnk9",
                 "Trpa1/Nos1/Efcab6",
                 "Selp/Fabp4/Vwf",
                 "Miat/Snhg11/Meg3",
                 "Olfm3/Tmem233/P2ry1",
                 "Sprr1a/Ecel1/Nts",
                 "Slc18a3/Pappa2/Runx3",
                 "Podxl/Depp1/Ptprb",
                 "Lox/Pkib/Kcnk3",
                 "Slc17a7/Hapln4/Lgi3",
                 "Trappc3l/Nptx1/Tafa1",
                 "Cd52/Lsp1/Rac2",
                 "Bmp3/Gal/Sost",
                 "Wfdc2/C1ql4/Rarres1",
                 "Gabra1/Chodl/Gabrb2",
                 "Thsd7b/Cacng5/Gata3",
                 "Hbb-bs/Bnc2/Ebf2",
                 "C1qc/Tyrobp/Lyz2",
                 "Gabra1/Gabrb2/Chrm2",
                 "Trpm8/Foxp2/Gsg1l",
                 "Smoc2/Ccl11/Lum",
                 "Ifitm1/Kcnj8/Rgs5",
                 "Slc6a2/Npy/Hand2",
                 "Amigo2/Glp1r/Npy2r",
                 "Olfr78/Runx3/Avpr1a")
seurat$genemarkers <- genemarkers[as.integer(seurat$cluster) + 1]
write.table(seurat, file = "leftright_deg_updated.txt", sep = "\t", quote = F, row.names = F)

