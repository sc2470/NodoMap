#Find differential expressed genes comparing fed and fasted

library(Seurat)

data<-readRDS("diet_nodomap_cleaned.RDS")
DefaultAssay(data) = "originalexp"
data <- SetIdent(data, value = "seurat_clusters")
levels(Idents(data))

#Subset sn-RNA seq inhouse data
table(paste0(data$Dataset, data$Nutr.cond))

inhouse <- subset(data, subset = Dataset == 'inhouse')
DefaultAssay(inhouse) = "originalexp"
levels(Idents(data))

#number of fedfasted cells in each clusters
table(paste0(inhouse@active.ident, inhouse$Nutr.cond))

##some groups have fewer than 3 cells
#7, 8, 18, 19, 20, 30, 31, 33, 34, 39, 42, 45, 46, 48, 49, 50, 51

#Fed and Fasted markers within each clusters
fedfastmarkers <- list()
for (i in c(0:6, 9:17, 21:29, 32, 35:38, 40:41, 43:44, 47)){
  print(i);fedfastmarkers[[i+1]] <- 
    FindMarkers(inhouse, ident.1 = "fast", ident.2 = "adlib", group.by = "Nutr.cond", subset.ident = i)
}


#Add gene names to fed and fasted markers
mouse_gene_names<-read.delim("mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

for (i in c(0:6, 9:17, 21:29, 32, 35:38, 40:41, 43:44, 47)){
  print(i);fedfastmarkers[[i+1]]$genesymbol <- 
    rownames(mouse_gene_names[match(rownames(fedfastmarkers[[i+1]]), mouse_gene_names$Gene.stable.ID),])
}

#Save the DEGs of each cluster
for (i in c(0:6, 9:17, 21:29, 32, 35:38, 40:41, 43:44, 47)){
  write.table(fedfastmarkers[[i+1]], print(paste0("fedfasted", i, ".txt")), 
              sep='\t', quote=F, col.names = NA)
}




# Seurat_FedFast_DEGs_barplot -----------------------------------------------------
library(ggplot2)
library(reshape2)
library(dplyr)

fedfastmarkers <- list.files(path = ".", pattern = "^fedfasted\\d+.txt")
fedfastmarkers <- do.call(rbind, Map("cbind", lapply(fedfastmarkers, read.delim), cluster=fedfastmarkers))
fedfastmarkers$cluster <- gsub("\\D","",fedfastmarkers$cluster)
colnames(fedfastmarkers)[1] <- "gene"

fedfastmarkers <- subset(fedfastmarkers, subset = p_val < 0.05)

degs <- fedfastmarkers %>%
  group_by(cluster) %>%
  dplyr::count(cluster)
up <- fedfastmarkers %>%
  group_by(cluster) %>%
  summarise(upregulated = sum(avg_log2FC>0))
down <- fedfastmarkers %>%
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

temp <- melt(temp, "clustermarkers")
temp <- subset(temp, variable %in% c("Fast_Up","Fast_Down"))
temp$clustermarkers <- factor(temp$clustermarkers, levels = unique(temp$clustermarkers))

split.row <- gsub("\\d","",temp$clustermarkers)
temp <- split.data.frame(temp, f = factor(split.row))

temp1 <- within(temp, rm("NGN", "JGN"))
temp1 <- do.call(rbind,temp1)
temp1$clustermarkers <- factor(temp1$clustermarkers, levels = unique(temp1$clustermarkers))

pdf("figure_fedfast_DEGs.pdf", width = 12, height = 4)
ggplot(temp$NGN, aes(x = clustermarkers, y = value, fill = variable)) +
  geom_col(position ="stack") + theme_classic() + 
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
  geom_col(position ="stack") + theme_classic() + 
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
  coord_fixed(ratio = 0.005)

ggplot(temp1, aes(x = clustermarkers, y = value, fill = variable)) +
  geom_col(position ="stack") + theme_classic() + 
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




# -------------------------------------------------------------------------


#Find differential expressed genes comparing left and right -------------
table(paste0(data$Dataset, data$Position))

#Subset Left & Right RNA seq data in Buchanan and inhouse
position <- subset(data, subset = Position != "unassigned")
DefaultAssay(position) = "originalexp"
levels(Idents(position))
table(paste0(position@active.ident, position$Position))

#Left & Right markers in each cluster (Buchanan and inhouse)
positionmarkers <- list()

for (i in c(0:19, 21:52)) {
  print(i);positionmarkers[[i+1]] <- 
    FindMarkers(position, ident.1 = "right", ident.2 = "left", group.by = "Position", subset.ident = i)
}

for (i in c(0:19, 21:52)){
  print(i);positionmarkers[[i+1]]$genesymbol <- 
    rownames(mouse_gene_names[match(rownames(positionmarkers[[i+1]]), mouse_gene_names$Gene.stable.ID),])
}

#Save the DEGs of each cluster
for (i in c(0:19, 21:52)){
  write.table(positionmarkers[[i+1]], print(paste0("leftright", i, ".txt")), sep='\t', quote=F, col.names = NA)
}



# Seurat_LeftRight_DEGs_barplot -------------------------------------------
leftrightmarkers <- list.files(path = ".", pattern = "^leftright\\d+.txt")
leftrightmarkers <- do.call(rbind, Map("cbind", lapply(leftrightmarkers, read.delim), cluster=leftrightmarkers))
leftrightmarkers$cluster <- gsub("\\D","",leftrightmarkers$cluster)
colnames(leftrightmarkers)[1] <- "gene"

leftrightmarkers <- subset(leftrightmarkers, subset = p_val < 0.05)

degs <- leftrightmarkers %>%
  group_by(cluster) %>%
  dplyr::count(cluster)
up <- leftrightmarkers %>%
  group_by(cluster) %>%
  summarise(upregulated = sum(avg_log2FC>0))
down <- leftrightmarkers %>%
  group_by(cluster) %>%
  summarise(downregulated = sum(avg_log2FC<0))
degs$Right_Up <- up$upregulated[match(degs$cluster, up$cluster)]
degs$Right_Down <- down$downregulated[match(degs$cluster, down$cluster)]

degs$cluster <- as.numeric(degs$cluster)
degs <- degs %>%
  arrange(cluster)

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



pdf("figure_leftright_DEGs.pdf", width = 8, height = 9)
ggplot(temp$NGN, aes(x = value, y = clustermarkers, fill = variable)) +
  geom_col(position = "stack") + theme_classic() + 
  scale_fill_manual(values = c("Right_Up" = "#AA3377", "Right_Down" = "#4477AA"), ##scale_fill_manual() matches fill colors to the levels of the variable
                    labels = c("Enriched in Right", "Enriched in Left")) + ##levels(factor(temp$NGN$variable))
  xlab("Number of DEGs") +
  ylab(NULL) + guides(fill= guide_legend(byrow = T, title = NULL)) +
  theme(axis.text = element_text(size = 19, color = "black", family = "Helvetica"),
        axis.title = element_text(size = 19),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 19),
        legend.spacing.y = unit(10, 'pt')) +
  coord_fixed(ratio = 200)

ggplot(temp$JGN, aes(x = value, y = clustermarkers, fill = variable)) +
  geom_col(position = "stack") + theme_classic() + 
  scale_fill_manual(values = c("Right_Up" = "#AA3377", "Right_Down" = "#4477AA"),
                    labels = c("Enriched in Right", "Enriched in Left")) +
  xlab("Number of DEGs") +
  ylab(NULL) + guides(fill= guide_legend(byrow = T, title = NULL)) +
  theme(axis.text = element_text(size = 19, color = "black", family = "Helvetica"),
        axis.title = element_text(size = 19),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 19),
        legend.spacing.y = unit(10, 'pt')) +
  coord_fixed(ratio = 50)

ggplot(temp1, aes(x = value, y = clustermarkers, fill = variable)) +
  geom_col(position = "stack") + theme_classic() + 
  scale_fill_manual(values = c("Right_Up" = "#AA3377", "Right_Down" = "#4477AA"),
                    labels = c("Enriched in Right", "Enriched in Left")) +
  xlab("Number of DEGs") +
  ylab(NULL) + guides(fill= guide_legend(byrow = T, title = NULL)) +
  theme(axis.text = element_text(size = 19, color = "black", family = "Helvetica"),
        axis.title = element_text(size = 19),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 19),
        legend.spacing.y = unit(10, 'pt'))
dev.off()














