library(Seurat)
library(ggplot2)
data_dir<-'/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/'
source(paste0(data_dir, "utility_functions.R"))
source(paste0(data_dir, "plotting_functions.R"))

load_plot_params()
load_colors()
getLongPalette = colorRampPalette(long_palette_strong)
getOkabeItoPalette = colorRampPalette(short_palette)

gene_names1<-read.delim(paste0(data_dir, "mouse_genes_names_Ens100_rmdup.txt"),row.names=2)
gene_names<-read.delim(paste0(data_dir, "mouse_genes_names_Ens100_rmdup.txt"),row.names=3)
genes_to_ens<-function(gene_symbol){
  return(as.character(gene_names[gene_symbol,1]))
}

##### Figure 1 and corresponding supp figs #####
nodomap<-readRDS(paste0(data_dir, 'diet_nodomap_cleaned.RDS'))


##Plot and save the NT and NP UMAPs
#NB something wrong with plotting quality
plot<-DimPlot(nodomap, 
              group.by = 'NPs',
              raster=TRUE,
              raster.dpi = c(rasterize_px,rasterize_px),
              cols=getOkabeItoPalette(20),
              pt.size = seurat_pt_size,
              label = F,
              repel = TRUE) + 
  NoAxes() + 
  theme(text = element_text(size=text_size), title = element_blank())
ggsave(filename = paste0(data_dir,"allcells_NPs.pdf"),
       plot = plot, "pdf",dpi=1200,width=450,height = 300,units="mm")
ggsave(filename = paste0(data_dir,"allcells_NPs.png"),
       plot = plot, "png",dpi=1200,width=450,height = 300,units="mm")

Idents(nodomap)<-'cell_identity'
neurons <- subset(nodomap, idents = grep("neuron", Idents(nodomap), value = TRUE, ignore.case = TRUE))

plot<-DimPlot(neurons, 
              group.by = 'NPs',
              raster=TRUE,
              raster.dpi = c(rasterize_px,rasterize_px),
              cols=getOkabeItoPalette(20),
              pt.size = seurat_pt_size,
              label = F,
              repel = TRUE) + 
  NoAxes() + 
  NoLegend() +
  theme(text = element_text(size=text_size), title = element_blank())
ggsave(filename = paste0(data_dir,"neurons_NPs_nolegend.pdf"),
       plot = plot, "pdf",dpi=1200,width=300,height = 300,units="mm")
ggsave(filename = paste0(data_dir,"neurons_NPs_nolegend.png"),
       plot = plot, "png",dpi=1200,width=300,height = 300,units="mm")

plot<-DimPlot(nodomap, 
              group.by = 'NTs',
              raster=TRUE,
              raster.dpi = c(rasterize_px,rasterize_px),
              cols=getOkabeItoPalette(4),
              pt.size = seurat_pt_size,
              label = F,
              repel = TRUE) + 
  NoAxes() + 
  theme(text = element_text(size=text_size), title = element_blank())
ggsave(filename = paste0(data_dir,"allcells_NTs.pdf"),
       plot = plot, "pdf",dpi=1200,width=450,height = 300,units="mm")
ggsave(filename = paste0(data_dir,"allcells_NTs.png"),
       plot = plot, "png",dpi=1200,width=450,height = 300,units="mm")

plot<-DimPlot(neurons, 
              group.by = 'NTs',
              raster=TRUE,
              raster.dpi = c(rasterize_px,rasterize_px),
              cols=getOkabeItoPalette(4),
              pt.size = seurat_pt_size,
              label = F,
              repel = TRUE) + 
  NoAxes() + 
  NoLegend()+
  theme(text = element_text(size=text_size), title = element_blank())
ggsave(filename = paste0(data_dir,"neurons_NTs_nolegend.pdf"),
       plot = plot, "pdf",dpi=1200,width=300,height = 300,units="mm")
ggsave(filename = paste0(data_dir,"neurons_NTs_nolegend.png"),
       plot = plot, "png",dpi=1200,width=300,height = 300,units="mm")

#Vlnplots of the genes of interest: 
plot<-VlnPlot(neurons, genes_to_ens(c('Slc17a6', 'Slc17a7', 'Dbh', 'Ddc', 'Th')), group.by = 'cluster_number', stack = T, cols = getOkabeItoPalette(5)) + 
  NoLegend() + 
  theme(axis.title = element_blank(), text = element_text(size=text_size))
ggsave(filename = paste0(data_dir,"neurons_NTs_vlnplot.pdf"),
       plot = plot, "pdf",dpi=1200,width=200,height = 200,units="mm")


plot<-VlnPlot(neurons, 
              genes_to_ens(c('Slc17a6', 'Slc17a7', 'Dbh', 'Ddc', 'Th', 'Gpha2', 'Adcyap1', 'Gal', 'Cartpt', 'Bdnf', 'Vip', 'Uts2b', 'Calca', 'Calcb', 'Tac1', 'Igf1', 'Nts', 'Ucn', 'Agrp', 'Npy')), 
              group.by = 'cluster_number', stack = T, cols = getOkabeItoPalette(20)) + 
  NoLegend() + 
  theme(axis.title = element_blank(), text = element_text(size=text_size))
ggsave(filename = paste0(data_dir,"neurons_NPs_NTs_vlnplot.pdf"),
       plot = plot, "pdf",dpi=1200,width=600,height = 200,units="mm")



