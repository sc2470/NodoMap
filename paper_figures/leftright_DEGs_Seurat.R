#Find differentially expressed genes comparing left and right
library(Seurat)

data<-readRDS("/home/rstudio/2023_April_Integration/nodose_integrated_PC30-60_res0.6-2/nodose-integrated-pcs-50-resolution-1.4.RDS")
DefaultAssay(data) = "originalexp"
data <- SetIdent(data, value = "integrated_snn_res.1.4")

table(paste0(data@meta.data$Experiment, data@meta.data$Position))

#Subset Left & Right RNA seq data in Buchanan and inhouse
position <- subset(data, subset = Position != "unassigned")
DefaultAssay(position) = "originalexp"
position <- SetIdent(position, value = "integrated_snn_res.1.4")

table(paste0(position@active.ident, position@meta.data$Position))

levels(position)

#Left & Right markers in each cluster (Buchanan and inhouse)
positionmarkers <- list()

for (i in c(0:52)) {
  print(i);positionmarkers[[i+1]] <- 
    FindMarkers(position, ident.1 = "right", ident.2 = "left", group.by = "Position", subset.ident = i)
}

#Add gene names to left and right markers
mouse_gene_names<-read.delim("../../mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}

for (i in c(0:52)){
  print(i);positionmarkers[[i+1]]$genesymbol <- 
    rownames(mouse_gene_names[match(rownames(positionmarkers[[i+1]]), mouse_gene_names$Gene.stable.ID),])
}

#Save the DEGs of each cluster
for (i in c(0:52)){
  write.table(positionmarkers[[i+1]], print(paste0("leftright", i, ".txt")), sep='\t', quote=F, col.names = NA)
}

saveRDS(positionmarkers, "nodose_leftrightmarkers_230811.RDS")

















