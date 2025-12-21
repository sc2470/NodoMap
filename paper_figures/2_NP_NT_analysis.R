library(Seurat)
library(ggplot2)
library(viridisLite)
library(viridis)
library(scales)
source('/cephfs/scratch/bruening_scratch/gdowsett/projects/nodomap/RNAseqSTUtils.R')

data_dir <- '/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/'
nucseq <- readRDS(paste0(data_dir, 'nodomap_cleaned_20251025.RDS'))

####### Gene functions 
gene_names1<-read.delim(paste0(data_dir, "mouse_genes_names_Ens100_rmdup.txt"),row.names=2)
gene_names<-read.delim(paste0(data_dir, "mouse_genes_names_Ens100_rmdup.txt"),row.names=3)
genes_to_ens<-function(gene_symbol){
  return(as.character(gene_names[gene_symbol,1]))
}

########## Look at neurotransmitter expression across clusters ##########
NT_genes<-c('Slc17a6', 'Slc17a7', 'Slc17a8', 'Slc32a1', 'Gad1', 'Gad2', 'Slc6a5', 'Slc6a9', 'Slc18a3', 'Chat', 'Slc5a7', 'Slc6a3', 'Ddc', 'Th', 'Pnmt', 'Dbh', 'Tph2', 'Slc6a4', 'Hdc')
DefaultAssay(nucseq)<-'originalexp'

NT<-calc_percent_cluster(nucseq, genes = genes_to_ens(NT_genes), metadata_name = 'cluster_number', assay = 'originalexp')

nt_model<-NT
nt_model[nt_model<30] <- 0
nt_model[nt_model>=30]<-1
nt_model<-nt_model[,1:19]
NT<-calc_NTs(nt_model, mouse = TRUE)

nucseq$NTs <- NT$NTs[ Idents(nucseq) ]

########### Look at neuropeptide expression across clusters ########
#read in NP and hormone list from thesis work 
np_genes_final<-read.delim(paste0(data_dir, '../../projects/nodomap/NP_mouse_curated_final.txt'))
np_genes_final<-np_genes_final$Gene
np_genes_final<-np_genes_final[-28]
#add in Gh and Bdnf 
np_genes_final<-c(np_genes_final, 'Bdnf', 'Gh')

nps<-calc_percent_cluster(nucseq, genes = genes_to_ens(np_genes_final), metadata_name = 'cluster_number', assay = 'originalexp')
np_pct<-nps[1:97]

test<-nps
test$Cluster<-rownames(test)
test <- test %>%
  mutate_at(vars(98:194), list(log1p))

ggplot(data = test, aes(x = Bdnf, y = reorder(Cluster, Bdnf), fill = AvEx_Bdnf)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(option = "viridis") +
  labs(x = "Percentage of Cells Expressing Pomc", y = "Cluster") +
  theme_minimal() +
  theme(axis.title = element_text(size = 18), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12))


np_pct$NPs <- apply(np_pct, 1, function(row) {
  np_genes <- colnames(np_pct)[row > 30 & !colnames(np_pct) %in% c("Adcyap1", 
                                                                   "Calca",
                                                                   "Calcb",
                                                                   "Gal", 
                                                                   "Tac1", 
                                                                   "Cartpt", 
                                                                   "Bdnf")]
  if (row['Adcyap1'] > 50) {
    np_genes <- c(np_genes, 'Adcyap1')
  }
  if (row['Calca'] > 50) {
    np_genes <- c(np_genes, 'Calca')
  }
  if (row['Calcb'] > 50) {
    np_genes <- c(np_genes, 'Calcb')
  }
  if (row['Gal'] > 40) {
    np_genes <- c(np_genes, 'Gal')
  }
  if (row['Tac1'] > 50) {
    np_genes <- c(np_genes, 'Tac1')
  }
  if (row['Cartpt'] > 60) {
    np_genes <- c(np_genes, 'Cartpt')
  }
  if (row['Bdnf'] > 50) {
    np_genes <- c(np_genes, 'Bdnf')
  }
  if (length(np_genes) > 0) {
    if (!is.na(row['NPs']) && row['NPs'] != '') {
      return(paste(row['NPs'], paste(np_genes, collapse = '|'), sep = '|'))
    } else {
      return(paste(np_genes, collapse = '|'))
    }
  } else {
    return(row['NPs'])
  }
})

nucseq$NPs <- np_pct$NPs[ Idents(nucseq) ]

saveRDS(nucseq, file.path(data_dir, "nodomap_cleaned_20251025.RDS"))