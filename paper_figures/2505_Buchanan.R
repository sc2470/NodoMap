###Compare NodoMap with Buchanan dataset
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

#Load data
data <- readRDS("diet_nodomap.RDS")

#Read Buchanan annotation
buchanan <- read.csv("GSE185173_nodoseMerged.meta 2.csv")

#genes_to_ens function
mouse_gene_names<-read.delim("../../mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

#extract metadata and cell names
x <- data@meta.data
x$original_barcodes<-rownames(x)
x_b <- x[x$Dataset == 'Buchanan',]

#update the rownames of x_b
#-1_2 = L_
#-1_3 = R_

x_b$cell_names <- rownames(x_b)

x_b <- x_b %>%
  mutate(cell_names = case_when(
    str_ends(cell_names, "-1_2") ~ str_c("L_", str_sub(cell_names, end = -5)),
    str_ends(cell_names, "-1_3") ~ str_c("R_", str_sub(cell_names, end = -5)),
    TRUE ~ cell_names
  ))

rownames(buchanan) <- buchanan$X
colnames(buchanan)[1] <- "cell_names"

test <- merge(x_b, buchanan, by = 'cell_names', all.x = TRUE)
x_meta <- merge(x, test, by = 'original_barcodes', all.x = TRUE)
x_meta <- x_meta[,c(1,96)]

data <- AddMetaData(data, metadata = x_meta$predicted.id, col.name = "Buchanan_annotation")


####Sankey plot with Buchanan vs NodoMap
annotations <- merge(x, test, by = 'original_barcodes', all.x = TRUE)
annotations <- annotations[,c(1,47,96)]

#subset to Buchanan data only
pattern <- "-1_[23]$"

annotation_buchanan <- annotations %>%
  filter(grepl(pattern, original_barcodes))

#Sankey plot
library(networkD3)

##update dataset into the appropriate format
rownames(annotation_buchanan) <- annotation_buchanan$original_barcodes
annotation_buchanan <- annotation_buchanan[,-1]
annotation_buchanan_connection <- data.frame(table(annotation_buchanan))
annotation_buchanan_connection <- annotation_buchanan_connection[annotation_buchanan_connection$Freq>0,]

##create a dataframe to list all cluster names and assigns a number to each
nodes <- data.frame(name=c(as.character(as.character(annotation_buchanan_connection$cluster_names.x)),
                           as.character(annotation_buchanan_connection$predicted.id)) %>% unique())
##add source and target column to the connection table which has the asssinged cluster numbers 
annotation_buchanan_connection$IDsource <- match(annotation_buchanan_connection$cluster_names.x,
                                                 nodes$name)-1
annotation_buchanan_connection$IDtarget <- match(annotation_buchanan_connection$predicted.id,
                                                 nodes$name)-1
#plot the sankey plot
sn <- sankeyNetwork(Links = annotation_buchanan_connection, #connection table
                    Nodes = nodes, #dataframe with the list of nodes
                    Source = "IDtarget", #column with the source numbers
                    Target = "IDsource", #column with the target numbers
                    Value = "Freq", #column with the sankey weights
                    NodeID = "name", #column with the source and target names
                    sinksRight=FALSE, #plot positioning
                    fontSize = 8, #node font size
                    fontFamily = 'Helvetica' #node font type
                    )

saveNetwork(sn, "../Scripts/2505_Buchanan_sankey_plot.html", selfcontained = TRUE)







