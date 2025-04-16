#After quality control
library(Seurat)

data <- readRDS("/home/rstudio/shared_data/nodose/latest-integration-230403/nodose-integrated.RDS")

DefaultAssay(data) <- 'integrated'

data <- RunPCA(data, npcs = 100, verbose = F)

DIMS<-c(30, 40, 50, 60)
resolution.list<-seq(0.6, 2.0, 0.2)

for (i in DIMS) {
  #data <- RunUMAP(data, reduction = "pca", dims = 1:i, umap.method = "UWOT", metric = "cosine")
  data <- FindNeighbors(data,reduction = "pca",dims = 1:i)
  for (j in resolution.list) {
    data<-FindClusters(data, resolution = j)
    print(paste0(i, ",", j))
    print(table(data$seurat_clusters))
    saveRDS(data, paste0("nodose-integrated-pcs-",i,"-resolution-",j,".RDS"))
  }
}

