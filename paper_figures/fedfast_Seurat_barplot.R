library(ggplot2)
library(reshape2)
library(dplyr)

seurat <- list.files(path = ".", pattern = "^fedfasted\\d+.txt")
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
degs$Fast_Up <- up$upregulated[match(degs$cluster, up$cluster)]
degs$Fast_Down <- down$downregulated[match(degs$cluster, down$cluster)]

degs$cluster <- as.numeric(degs$cluster)
degs <- degs %>%
  arrange(cluster)

temp <- degs[,c(1,3,4)]
temp2 <- data.frame(cluster = c(7,8,18,19,30,31,33,34,39,42,45,46,48,49,50,51,52), Fast_Up = 0, Fast_Down = 0)
temp <- rbind(as.data.frame(temp), temp2)
temp<-temp[order(temp$cluster),]

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

temp$clustermarkers <- clustermarkers
temp$Fast_Down <- temp$Fast_Down*-1

temp <- melt(temp, "clustermarkers")
temp <- subset(temp, variable %in% c("Fast_Up","Fast_Down"))
temp$clustermarkers <- factor(temp$clustermarkers, levels = unique(temp$clustermarkers))

split.row <- gsub("\\d","",temp$clustermarkers)
temp <- split.data.frame(temp, f = factor(split.row))

temp1 <- within(temp, rm("NGN", "JGN"))
temp1 <- do.call(rbind,temp1)
temp1$clustermarkers <- factor(temp1$clustermarkers, levels = unique(temp1$clustermarkers))

pdf("figure_Seurat_DEGs.pdf", width = 12, height = 4)
ggplot(temp$NGN, aes(x = clustermarkers, y = value, fill = variable)) +
  geom_col(stat="stack") + theme_classic() + 
  scale_fill_manual(values = c("Fast_Up" = "#BB5566", "Fast_Down" = "#004488"), 
                    labels=c("Upregulated in Fasting", "Downregulated in Fasting")) +
  xlab(NULL) +
  ylab("Number of DEGs") + guides(fill= guide_legend(byrow = T, title = NULL)) +
  theme(axis.text = element_text(size = 21, color = "black", family = "Helvetica"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 21),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 21),
        legend.spacing.y = unit(10, 'pt'))


ggplot(temp$JGN, aes(x = clustermarkers, y = value, fill = variable)) +
  geom_col(stat="stack", width = 0.9) + theme_classic() + 
  scale_fill_manual(values = c("Fast_Up" = "#BB5566", "Fast_Down" = "#004488"),
                    labels=c("Upregulated in Fasting", "Downregulated in Fasting")) +
  xlab(NULL) +
  ylab("Number of DEGs") + guides(fill= guide_legend(byrow = T, title = NULL)) +
  theme(axis.text = element_text(size = 21, color = "black", family = "Helvetica"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 21),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 21),
        legend.spacing.y = unit(10, 'pt')) +
  coord_fixed(ratio = 0.007)

ggplot(temp1, aes(x = clustermarkers, y = value, fill = variable)) +
  geom_col(stat="stack") + theme_classic() + 
  scale_fill_manual(values = c("Fast_Up" = "#BB5566", "Fast_Down" = "#004488"),
                    labels=c("Upregulated in Fasting", "Downregulated in Fasting")) +
  xlab(NULL) +
  ylab("Number of DEGs") + guides(fill= guide_legend(byrow = T, title = NULL)) +
  theme(axis.text = element_text(size = 21, color = "black", family = "Helvetica"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 21),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 21),
        legend.spacing.y = unit(10, 'pt'))
dev.off()








