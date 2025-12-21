library(ggplot2)
library(reshape2)
library(dplyr)
library(scCustomize)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(stringr)
library(pheatmap)

# Extract the neuronal cluster annotations------------------------------------------
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


sodium_channel <- c("Nav1.8",
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


fibre_type <- c("Unmyelinated",
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


sensor_type <- c("Mechanosensor/Nocisensor",
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


organ_innervated <- c("Broad projection",
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

dt <- data.frame(Clusters = neurons,
                 Sodium_channel = sodium_channel,
                 Fibre_type = fibre_type,
                 Sensor_type = sensor_type,
                 Organ_projection = organ_innervated)

# Extract the DEGs info ---------------------------------------------------
###Overnight fasting vs Ad libitum
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
                    "SGC8",
                    "JGN1",
                    "NGN6",
                    "GC2",
                    "JGN2",
                    "FB3",
                    "NGN7",
                    "NGN8",
                    "EC2",
                    "NGN9",
                    "NGN10",
                    "NGN11",
                    "NGN12",
                    "EC3",
                    "NGN13",
                    "NGN14",
                    "JGN3",
                    "HC2",
                    "NGN15",
                    "JGN4",
                    "NGN16",
                    "NGN17",
                    "FB4",
                    "GC3",
                    "NGN18",
                    "JGN5",
                    "MGC5",
                    "MGC6",
                    "NGN19",
                    "NGN20",
                    "NGN21")

temp$clustermarkers <- clustermarkers
temp$Fast_Down <- temp$Fast_Down*-1
temp<-temp[order(temp$clustermarkers),]
temp <- temp[c(13:17,24:45),]

#Tag fasting DEGs to annotations
dt$Upregulated_in_Fasting <- temp$Fast_Up[match(dt$Clusters,temp$clustermarkers)]
dt$Downregulated_in_Fasting <- temp$Fast_Down[match(dt$Clusters,temp$clustermarkers)]

###Left and right nodose ganglia 
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

degs$clustermarkers <- clustermarkers
degs$Right_Down <- degs$Right_Down*-1
degs <- degs[order(degs$clustermarkers),]
degs <- degs[c(13:17,24:45),]

#Tag position DEGs to annotations
dt$Enriched_in_Right <- degs$Right_Up[match(dt$Clusters,degs$clustermarkers)]
dt$Enriched_in_Left <- degs$Right_Down[match(dt$Clusters,degs$clustermarkers)]

#Heatmap
rownames(dt) <- dt$Clusters
dt$Clusters <- NULL

degs <- dt[,c(5:8)]
classification <- dt[,c(1:4)]

col_class <- list(Sodium_channel = brewer.pal(8, "Set2")[2:4],
                  Fibre_type = brewer.pal(8, "Set1")[1:3],
                  Sensor_type = brewer.pal(8, "Dark2")[1:3],
                  Organ_projection = brewer.pal(8, "Accent")[1:7])

names(col_class$Sodium_channel) <- unique(classification$Sodium_channel)
names(col_class$Fibre_type) <- unique(classification$Fibre_type)
names(col_class$Sensor_type) <- unique(classification$Sensor_type)
names(col_class$Organ_projection) <- unique(classification$Organ_projection)

pdf("deg_anno_heatmap.pdf", width = 10, height = 15)
pheatmap(degs, color = colorRampPalette(rev(brewer.pal(n=9, name = "RdYlBu")))(100),
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




