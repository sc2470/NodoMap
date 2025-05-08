#This script is to compare our clustering and nucseq annotation to publicly available datasets
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

#1. read in data
data<-readRDS('~/OneDrive - University of Cambridge/NodoMap/datasets/diet_nodomap.RDS')

#Read in Kupari annotations
kupari<-read.delim('~/OneDrive - University of Cambridge/NodoMap/cluster-annotations-lit-comparisons/Kupari_GSE124312_Vagus_AnnotationTable.txt')

#genes_to_ens function
gene_names1<-read.delim("~/Dropbox/Thesis/Chapter2/Hindbrain/Reference_tables/mouse_genes_names_Ens100_rmdup.txt",row.names=2)
gene_names<-read.delim("~/Dropbox/Thesis/Chapter2/Hindbrain/Reference_tables/mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol){
  return(as.character(gene_names[gene_symbol,1]))
}

#extract metadata and cell names from datasets 
x<-data@meta.data
x$original_barcodes<-rownames(x)
x_k <- x[x$Dataset == 'Kupari',]


#For x_k, change the rownames accordingly: 
#-1_5 = 003_
#-1_6 = 005_
#-1_7 = 102_
#-1_8 = 103_

x_k$cell_names<-rownames(x_k)

x_k <- x_k %>%
  mutate(cell_names = case_when(
    str_ends(cell_names, "-1_5") ~ str_c("003_", str_sub(cell_names, end = -5)),
    str_ends(cell_names, "-1_6") ~ str_c("005_", str_sub(cell_names, end = -5)),
    str_ends(cell_names, "-1_7") ~ str_c("102_", str_sub(cell_names, end = -5)),
    str_ends(cell_names, "-1_8") ~ str_c("103_", str_sub(cell_names, end = -5)),
    TRUE ~ cell_names
  ))

kupari$cell_names <- rownames(kupari)

test<-merge(x_k, kupari, by = 'cell_names', all.x = TRUE)
x_meta<-merge(x, test, by = 'original_barcodes', all.x = TRUE)
x_meta<-x_meta[,c(1,99)]

data<-AddMetaData(data, metadata = x_meta$Cluster, col.name = 'Kupari_annotation')


######## Sankey plot of Kupari vs Our annotations, in just the Kupari cells and in all cells 

annotations<-merge(x, test, by = 'original_barcodes', all.x = TRUE)
annotations<-annotations[,c(1,47,99)]

#subset to kupari cells only 
pattern <- "-1_[5678]$"

annotation_kupari <- annotations %>% 
  filter(grepl(pattern, original_barcodes))

#Lets learn how to make Sankey plots 

#Package used
library(networkD3)

#first, we need to change the dataset into the appropriate format
rownames(annotation_kupari)<-annotation_kupari$original_barcodes
annotation_kupari<-annotation_kupari[,-1]
annotation_kupari_connection<-data.frame(table(annotation_kupari))
annotation_kupari_connection<-annotation_kupari_connection[annotation_kupari_connection$Freq>0,]

#Create a dataframe which lists all cluster names and assigns a number to each
nodes <- data.frame(name=c(as.character(annotation_kupari_connection$cluster_names.x),
                           as.character(annotation_kupari_connection$Cluster)) %>% unique())
#add source and target columns to the connection table which has the asssinged cluster numbers 
annotation_kupari_connection$IDsource <- match(annotation_kupari_connection$cluster_names.x, nodes$name)-1
annotation_kupari_connection$IDtarget <- match(annotation_kupari_connection$Cluster, nodes$name)-1

#plot the sankey plot
sn<-sankeyNetwork(Links = annotation_kupari_connection, #connection table
              Nodes = nodes, #data frame with list of nodes 
              Source = "IDtarget", # column with the source numbers 
              Target = "IDsource", #column with the target numbers 
              Value = "Freq",  #column with the sankey weights 
              NodeID = "name", #column with the source and target names 
              sinksRight=FALSE, #plot positioning
              fontSize = 8, #node font size
              fontFamily = 'Helvetica' #node font type
              )
saveNetwork(sn, '~/OneDrive - University of Cambridge/NodoMap/Scripts/sankey-plot.html', selfcontained = TRUE)


#Look to see the porportion of Kupari cells in each of our own clusters 

kupari_proportion<-data.frame(table(subset(data, Dataset == 'Kupari')$cluster_names))
all_clusters<-data.frame(table(data$cluster_names))
kupari_proportion<-merge(kupari_proportion, all_clusters, by = 'Var1')
kupari_proportion$proportion<-kupari_proportion$Freq.x/kupari_proportion$Freq.y*100

## looking at the cells which have been put into clusters which are different cell types between the two clustering strategies 
subset_connection <- annotation_kupari_connection[annotation_kupari_connection$Freq %in% c(1, 2), ]
misc <- merge(annotation_kupari, subset_connection,
              by.x = c("cluster_names.x", "Cluster"),
              by.y = c("cluster_names.x", "Cluster"))
## we have 55 cells that are mapped between clusters with a freq of 1 or 2. 
## we look at this list to identify first, cells that have different cell type identity between the two groups 
misc_cellchange<-misc[c(1:2, 5:7, 9, 11, 18, 28),]
## total of 9 cells that have different cell identity across the two clusterings. 
## subset and look at cell type marker genes for these cells 
cells<-as.character(misc_cellchange$cellID)

doublets<-subset(data, cells = cells)

#Look at cell type marker genes 
Idents(data)<-'cell_identity'
celltype_markers<-FindAllMarkers(data, only.pos = T)
celltype_markers$genename <- gene_names1$Gene.name[match(celltype_markers$gene, sub("\\..*", "", gene_names1$Gene.stable.ID))]
celltype_markers$specificity_score <- ((celltype_markers$pct.1 + 0.01) / (celltype_markers$pct.2 + 0.01)) * celltype_markers$avg_log2FC

neuron<-subset(data, idents = 'Neuron')
Idents(neuron)<-'cell_types'
neuron@meta.data$neuron_type <- substr(neuron@meta.data$cell_types, 1, 3)
neuron_markers<-FindAllMarkers(neuron, only.pos = T)
neuron_markers$genename <- gene_names1$Gene.name[match(neuron_markers$gene, sub("\\..*", "", gene_names1$Gene.stable.ID))]
neuron_markers$specificity_score <- ((neuron_markers$pct.1 + 0.01) / (neuron_markers$pct.2 + 0.01)) * neuron_markers$avg_log2FC

celltype_genes<-FetchData(doublets, vars = genes_to_ens(c('Phox2b', 'P2rx2', 'Hoxb5', 
                                                          'Prdm12', 'Tmem45b', 'Prrxl1', 
                                                          'Ly6c1', 'Cldn5', 'Ly6a', 'Emcn', 
                                                          'Tyrobp', 'C1qc', 'C1qa', 'Ptprc', 
                                                          'Smoc2', 'Ccl11', 'Dcn', 'Ebf2')))
colnames(celltype_genes)<-c('Phox2b', 'P2rx2', 'Hoxb5', 
                            'Prdm12', 'Tmem45b', 'Prrxl1', 
                            'Ly6c1', 'Cldn5', 'Ly6a', 'Emcn', 
                            'Tyrobp', 'C1qc', 'C1qa', 'Ptprc', 
                            'Smoc2', 'Ccl11', 'Dcn', 'Ebf2'
)
celltype_genes$cellID<-rownames(celltype_genes)

misc_cellchange<-merge(misc_cellchange, celltype_genes, by = 'cellID')
write.table(misc_cellchange, '~/OneDrive - University of Cambridge/NodoMap/cluster-annotations-lit-comparisons/kupari_misc_cells_gene_exp.txt', sep = '\t', row.names = F)

