#After run multi-PC/resolution
#Calculate average silhouette width (ASW) of each cluster under different PC&Res
#ASW indicates the within-cluster distance between the cells, 
#and the between-cluster distance of the cells to the closetest clusters

#Randomly select 25000 cells at pc30-60 and res0.6-2 and calculate their sil_width for 5 times
#Calculate the average silhouette width of each run

library(Seurat)
library(cluster)

setwd("/home/rstudio/2023_April_Integration/nodose_integrated_PC30-60_res0.6-2")
res.list<-seq(0.6, 2.0, 0.2)
sil_results<-vector("numeric",length=5)

#PC30

for (i in res.list) {
  data <- readRDS(paste0("nodose-integrated-pcs-30-resolution-",i,".RDS"))
  for (j in 1:5) {
    print(paste0("j = ",j))
    set.seed(j)
    print("subsetting")
    cellstosubset<-(sample(colnames(data),25000))
    #print(cellstosubset)
    nodose_sub<-subset(data,cells = cellstosubset)
    print("calculating dist matrix")
    dist.matrix<-dist(x=Embeddings(nodose_sub[["pca"]])[,1:50])
    print(head(dist.matrix))
    clusters<-nodose_sub@active.ident
    print("calculating silhouette width")
    sil<-silhouette(as.numeric(as.factor(clusters)),dist=dist.matrix)
    print(head(sil))
    print("clearing memory")
    gc()
    print("asw")
    print(summary(sil[,3]))
    sil_results[j]<-mean(sil[,3])
    write.table(sil, paste0("nodose-integrated-pcs-30-resolution-",i,"run",j,".sw.txt"),sep="\t",quote=F,col.names=NA)
  }
  write.table(sil_results,paste0("nodose-integrated-pcs-30-resolution-",i,".asw.txt"),sep="\t",quote=F,col.names=NA)
}



#PC40

for (i in res.list) {
  data <- readRDS(paste0("nodose-integrated-pcs-40-resolution-",i,".RDS"))
  for (j in 1:5) {
    print(paste0("j = ",j))
    set.seed(j)
    print("subsetting")
    cellstosubset<-(sample(colnames(data),25000))
    #print(cellstosubset)
    nodose_sub<-subset(data,cells = cellstosubset)
    print("calculating dist matrix")
    dist.matrix<-dist(x=Embeddings(nodose_sub[["pca"]])[,1:50])
    print(head(dist.matrix))
    clusters<-nodose_sub@active.ident
    print("calculating silhouette width")
    sil<-silhouette(as.numeric(as.factor(clusters)),dist=dist.matrix)
    print(head(sil))
    print("clearing memory")
    gc()
    print("asw")
    print(summary(sil[,3]))
    sil_results[j]<-mean(sil[,3])
    write.table(sil, paste0("nodose-integrated-pcs-40-resolution-",i,"run",j,".sw.txt"),sep="\t",quote=F,col.names=NA)
  }
  write.table(sil_results,paste0("nodose-integrated-pcs-40-resolution-",i,".asw.txt"),sep="\t",quote=F,col.names=NA)
}

#PC50

for (i in res.list) {
  data <- readRDS(paste0("nodose-integrated-pcs-50-resolution-",i,".RDS"))
  for (j in 1:5) {
    print(paste0("j = ",j))
    set.seed(j)
    print("subsetting")
    cellstosubset<-(sample(colnames(data),25000))
    #print(cellstosubset)
    nodose_sub<-subset(data,cells = cellstosubset)
    print("calculating dist matrix")
    dist.matrix<-dist(x=Embeddings(nodose_sub[["pca"]])[,1:50])
    print(head(dist.matrix))
    clusters<-nodose_sub@active.ident
    print("calculating silhouette width")
    sil<-silhouette(as.numeric(as.factor(clusters)),dist=dist.matrix)
    print(head(sil))
    print("clearing memory")
    gc()
    print("asw")
    print(summary(sil[,3]))
    sil_results[j]<-mean(sil[,3])
    write.table(sil, paste0("nodose-integrated-pcs-50-resolution-",i,"run",j,".sw.txt"),sep="\t",quote=F,col.names=NA)
  }
  write.table(sil_results,paste0("nodose-integrated-pcs-50-resolution-",i,".asw.txt"),sep="\t",quote=F,col.names=NA)
}

#PC60

for (i in res.list) {
  data <- readRDS(paste0("nodose-integrated-pcs-60-resolution-",i,".RDS"))
  for (j in 1:5) {
    print(paste0("j = ",j))
    set.seed(j)
    print("subsetting")
    cellstosubset<-(sample(colnames(data),25000))
    #print(cellstosubset)
    nodose_sub<-subset(data,cells = cellstosubset)
    print("calculating dist matrix")
    dist.matrix<-dist(x=Embeddings(nodose_sub[["pca"]])[,1:50])
    print(head(dist.matrix))
    clusters<-nodose_sub@active.ident
    print("calculating silhouette width")
    sil<-silhouette(as.numeric(as.factor(clusters)),dist=dist.matrix)
    print(head(sil))
    print("clearing memory")
    gc()
    print("asw")
    print(summary(sil[,3]))
    sil_results[j]<-mean(sil[,3])
    write.table(sil, paste0("nodose-integrated-pcs-60-resolution-",i,"run",j,".sw.txt"),sep="\t",quote=F,col.names=NA)
  }
  write.table(sil_results,paste0("nodose-integrated-pcs-60-resolution-",i,".asw.txt"),sep="\t",quote=F,col.names=NA)
}


#Calculate the average silhouette width of each cluster
data <- list.files(path = ".", pattern = "[1-9].sw.txt")

myaggregate<-function(data)
{
  aggregate(data$sil_width, list(data$cluster), FUN=mean)
}


lapply(data,FUN = function(x){
  y<-read.delim(x,sep="\t",row.names = 1)
  z<-myaggregate(y)
  write.table(z,paste0(x,"_pergroup_asw.txt"),sep="\t",quote=F,col.names=NA)
})


















