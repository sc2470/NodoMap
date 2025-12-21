### Inspect nodomap output 
library(Seurat)

data_dir <- '/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/'
data <-readRDS(paste0(data_dir, 'nodomap_cleaned_20251025.RDS'))

data <- DietSeurat(
  data,
  assays = 'originalexp',
  dimreducs = 'umap',
  graphs = c('integrated_nn', 'integrated_snn'),
  misc = FALSE
)
DimPlot(data)

saveRDS(data, paste0(data_dir, 'diet_nodomap_cleaned.RDS'))
