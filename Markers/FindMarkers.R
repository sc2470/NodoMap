#After identified the most fitted PC and resolution
library(Seurat)
library(dplyr)

data <- readRDS("/home/rstudio/2023_April_Integration/nodose_integrated_PC30-60_res0.6-2/nodose-integrated-pcs-50-resolution-1.4.RDS")
DefaultAssay(data) = "originalexp"
data <- SetIdent(data, value = "integrated_snn_res.1.4")

#Find markers of every clusters, only show positive ones
data.marker <- FindAllMarkers(data, only.pos = T, min.pct = 0.25)

View(data.marker)

mouse_gene_names<-read.delim("../mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol)+   
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

data.marker$genesymbol <- rownames(mouse_gene_names[match(data.marker$gene, mouse_gene_names$Gene.stable.ID),])

write.table(data.marker, file = "pc50-resolution-1.4-markers-230613.txt", sep = "\t", quote = F, col.names = NA)

data.marker %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)

#Label clusters with marker genes
markers <- read.delim("/home/rstudio/2023_April_Integration/pc50-resolution-1.4-markers-230613.txt")

markers$prop <- markers$avg_log2FC*(markers$pct.1 + 0.01)/(markers$pct.2 + 0.01)

top3 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = prop)

