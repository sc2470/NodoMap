#Find differentially expressed genes comparing fed and fasted

library(Seurat)

data<-readRDS("/home/rstudio/2023_April_Integration/nodose_integrated_PC30-60_res0.6-2/nodose-integrated-pcs-50-resolution-1.4.RDS")
DefaultAssay(data) = "originalexp"
data <- SetIdent(data, value = "integrated_snn_res.1.4")

table(paste0(data@active.ident, data@meta.data$Nutr.cond))

#Subset sn-RNA seq inhouse data
table(data@meta.data$Dataset)

inhouse <- subset(data, subset = Dataset == 'inhouse')
DefaultAssay(inhouse) = "originalexp"
inhouse <- SetIdent(inhouse, value = "integrated_snn_res.1.4")

#number of fedfasted cells in each clusters
table(paste0(inhouse@active.ident, inhouse@meta.data$Nutr.cond))

##some groups have fewer than 3 cells
#7, 8, 18, 19, 30, 31, 33, 34, 39, 42, 45, 46, 48, 49, 50, 52 

#Fed and Fasted markers within each clusters
fedfastmarkers <- list()
for (i in c(0:6, 9:17, 20:29, 32, 35:38, 40:41, 43:44, 47)){
  print(i);fedfastmarkers[[i+1]] <- 
    FindMarkers(inhouse, ident.1 = "fast", ident.2 = "adlib", group.by = "Nutr.cond", subset.ident = i)
}


#Add gene names to fed and fasted markers
mouse_gene_names<-read.delim("../../mouse_genes_names_Ens100_rmdup.txt",row.names=3)
genes_to_ens<-function(gene_symbol) + 
  { return(as.character(mouse_gene_names[gene_symbol,1]))}


for (i in c(0:6, 9:17, 20:29, 32, 35:38, 40:41, 43:44, 47)){
  print(i);fedfastmarkers[[i+1]]$genesymbol <- 
    rownames(mouse_gene_names[match(rownames(fedfastmarkers[[i+1]]), mouse_gene_names$Gene.stable.ID),])
}


#Save the fed and fasted markers of each clusters
#for (i in c(0:6, 9:17, 20:29, 32, 35:38, 40:41, 43:44, 47)){
#  fedfastmarkers[[i+1]] <- subset(fedfastmarkers[[i+1]], subset=fedfastmarkers[[i+1]]$p_val < 0.05)
#  write.table(fedfastmarkers[[i+1]], print(paste0("fedfasted", i, ".txt")), sep='\t', quote=F, col.names = NA)
#}

saveRDS(fedfastmarkers, "nodose_fedfastmarkers_230808.RDS")

#Save for pathway analysis
for (i in c(0:6, 9:17, 20:29, 32, 35:38, 40:41, 43:44, 47)){
  write.table(fedfastmarkers[[i+1]], print(paste0("fedfasted", i, ".txt")), sep='\t', quote=F, col.names = NA)
}




