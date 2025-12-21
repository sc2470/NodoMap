
# ###Compare NodoMap with Kupari dataset ----------------------------------
###Compare NodoMap with Buchanan dataset
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

#Load data
data <- readRDS("../data_objects/diet_nodomap_metadata_updated_2510.rds")


#Read in Kupari annotations
kupari<-read.delim('../../cluster-annotations-lit-comparisons/Kupari_GSE124312_Vagus_AnnotationTable.txt')

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
x_meta<-x_meta[,c(1,89)]



######## Sankey plot of Kupari vs Our annotations, in just the Kupari cells and in all cells 

annotations<-merge(x, test, by = 'original_barcodes', all.x = TRUE)
annotations<-annotations[,c(1,38,89)]

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
nodes <- data.frame(name=c(as.character(annotation_kupari_connection$cluster_markers.x),
                           as.character(annotation_kupari_connection$Cluster)) %>% unique())
#add source and target columns to the connection table which has the asssinged cluster numbers 
annotation_kupari_connection$IDsource <- match(annotation_kupari_connection$cluster_markers.x, nodes$name)-1
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
saveNetwork(sn, 'Kupari_sankey_plot.html', selfcontained = TRUE)












