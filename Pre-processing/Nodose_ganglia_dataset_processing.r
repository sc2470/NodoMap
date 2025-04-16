library(scater)
library(DropletUtils)
library(scDblFinder)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(hdf5r)
library(tidyr)

#bai
setwd(dir = "../../../genehunter/nodose")
bai<-read10xCounts("bai2019/filtered_feature_bc_matrix", col.names = T)
mt_genes<-rownames(bai[grep(pattern = "mt-", x = bai@rowRanges@elementMetadata$Symbol),])
mt_genes<-mt_genes[-1]
bai<-addPerCellQC(bai, subset = list(Mito = mt_genes))
plotColData(bai, x = "sum", y = "subsets_Mito_percent")
x<-plotColData(bai, y = "sum")
y<-plotColData(bai, y = "detected")
z<-plotColData(bai, y = "subsets_Mito_percent")
x + y + z
summary(bai$sum)
summary(bai$detected)
summary(bai$subsets_Mito_percent)
keep.total<-bai$sum >1000
keep.n<-bai$detected>200
keepmito<-bai$subsets_Mito_percent<10
bai<-bai[,keep.total & keep.n & keepmito]
set.seed(1)
bai<-scDblFinder(bai)
#log normalise to allow for conversion to seurat object 
bai<-logNormCounts(bai)
bai1_seurat<-as.Seurat(bai)
#add meta data
bai1_seurat<-AddMetaData(bai1_seurat, "Bai", col.name = "Dataset")
bai1_seurat<-AddMetaData(bai1_seurat, "Bai_10X", col.name = "Experiment")
bai1_seurat<-AddMetaData(bai1_seurat, "6 weeks", col.name = "Age")
bai1_seurat<-AddMetaData(bai1_seurat, "WT", col.name = "Strain")
# remove doublets
bai1_seurat<-subset(bai1_seurat, subset= scDblFinder.class == "singlet")
VlnPlot(bai1_seurat, features = c("nFeature_originalexp", "nCount_originalexp", "subsets_Mito_percent"), ncol = 3)
saveRDS(bai1_seurat, "../../mouse-nuc-seq/nodose/data/processed/bai-processed-220711.RDS")


#kupari 102 dataset
kupari102<-read10xCounts("kupari2019/10X_17_102/filtered_feature_bc_matrix", col.names = T)
mt_genes<-rownames(kupari102[grep(pattern = "mt-", x = kupari102@rowRanges@elementMetadata$Symbol),])
mt_genes<-mt_genes[-1]
kupari102<-addPerCellQC(kupari102, subset = list(Mito = mt_genes))
plotColData(kupari102, x = "sum", y = "subsets_Mito_percent")
x<-plotColData(kupari102, y = "sum")
y<-plotColData(kupari102, y = "detected")
z<-plotColData(kupari102, y = "subsets_Mito_percent")
x + y + z
summary(kupari102$sum)
summary(kupari102$detected)
summary(kupari102$subsets_Mito_percent)
keep.total<-kupari102$sum >1000
keep.top<-kupari102$sum < 100000
keep.n<-kupari102$detected>200
keepmito<-kupari102$subsets_Mito_percent<10
kupari102<-kupari102[,keep.total & keep.n & keepmito & keep.top]
set.seed(1)
kupari102<-scDblFinder(kupari102)
#log normalise to allow for conversion to seurat object 
kupari102<-logNormCounts(kupari102)
kupari102_seurat<-as.Seurat(kupari102)

#add meta data
kupari102_seurat<-AddMetaData(kupari102_seurat, "Kupari", col.name = "Dataset")
kupari102_seurat<-AddMetaData(kupari102_seurat, "Kupari_102", col.name = "Experiment")
kupari102_seurat<-AddMetaData(kupari102_seurat, "5 weeks", col.name = "Age")
kupari102_seurat<-AddMetaData(kupari102_seurat, "WT", col.name = "Strain")
# remove doublets
kupari102_seurat<-subset(kupari102_seurat, subset= scDblFinder.class == "singlet")
VlnPlot(kupari102_seurat, features = c("nFeature_originalexp", "nCount_originalexp", "subsets_Mito_percent"), ncol = 3)
saveRDS(kupari102_seurat, "../../mouse-nuc-seq/nodose/data/processed/kupari102-seurat-processed-220711.RDS")

#kupari 103 dataset
kupari103<-read10xCounts("kupari2019/10X_17_103/filtered_feature_bc_matrix", col.names = T)
mt_genes<-rownames(kupari103[grep(pattern = "mt-", x = kupari103@rowRanges@elementMetadata$Symbol),])
mt_genes<-mt_genes[-1]
kupari103<-addPerCellQC(kupari103, subset = list(Mito = mt_genes))
plotColData(kupari103, x = "sum", y = "subsets_Mito_percent")

x<-plotColData(kupari103, y = "sum")
y<-plotColData(kupari103, y = "detected")
z<-plotColData(kupari103, y = "subsets_Mito_percent")
x + y + z
summary(kupari103$sum)
summary(kupari103$detected)
summary(kupari103$subsets_Mito_percent)

keep.total<-kupari103$sum >1000
keep.top<-kupari103$sum < 75000
keep.n<-kupari103$detected>200
keepmito<-kupari103$subsets_Mito_percent<10

kupari103<-kupari103[,keep.total & keep.n & keepmito & keep.top]

set.seed(1)
kupari103<-scDblFinder(kupari103)
#log normalise to allow for conversion to seurat object 
kupari103<-logNormCounts(kupari103)
kupari103_seurat<-as.Seurat(kupari103)
#add meta data
kupari103_seurat<-AddMetaData(kupari103_seurat, "Kupari", col.name = "Dataset")
kupari103_seurat<-AddMetaData(kupari103_seurat, "Kupari_103", col.name = "Experiment")
kupari103_seurat<-AddMetaData(kupari103_seurat, "5 weeks", col.name = "Age")
kupari103_seurat<-AddMetaData(kupari103_seurat, "WT", col.name = "Strain")
# remove doublets
kupari103_seurat<-subset(kupari103_seurat, subset= scDblFinder.class == "singlet")

VlnPlot(kupari103_seurat, features = c("nFeature_originalexp", "nCount_originalexp", "subsets_Mito_percent"), ncol = 3)
saveRDS(kupari103_seurat, "../../mouse-nuc-seq/nodose/data/processed/kupari103-seurat-processed-220711.RDS")

#kupari 003 dataset
kupari003<-read10xCounts("kupari2019/10X_18_003/filtered_feature_bc_matrix", col.names = T)
mt_genes<-rownames(kupari003[grep(pattern = "mt-", x = kupari003@rowRanges@elementMetadata$Symbol),])
mt_genes<-mt_genes[-1]
kupari003<-addPerCellQC(kupari003, subset = list(Mito = mt_genes))
plotColData(kupari003, x = "sum", y = "subsets_Mito_percent")
#plotting QC data
x<-plotColData(kupari003, y = "sum")
y<-plotColData(kupari003, y = "detected")
z<-plotColData(kupari003, y = "subsets_Mito_percent")
x + y + z
summary(kupari003$sum)
summary(kupari003$detected)
summary(kupari003$subsets_Mito_percent)
#removing low quality cells 
keep.total<-kupari003$sum >1000
keep.top<-kupari003$sum < 75000
keep.n<-kupari003$detected>400
keepmito<-kupari003$subsets_Mito_percent<10
kupari003<-kupari003[,keep.total & keep.n & keepmito & keep.top]
set.seed(1)
kupari003<-scDblFinder(kupari003)
#log normalise to allow for conversion to seurat object 
kupari003<-logNormCounts(kupari003)
kupari003_seurat<-as.Seurat(kupari003)
#add meta data
kupari003_seurat<-AddMetaData(kupari003_seurat, "Kupari", col.name = "Dataset")
kupari003_seurat<-AddMetaData(kupari003_seurat, "Kupari_003", col.name = "Experiment")
kupari003_seurat<-AddMetaData(kupari003_seurat, "5 weeks", col.name = "Age")
kupari003_seurat<-AddMetaData(kupari003_seurat, "Vglut2Cre-Tomato", col.name = "Strain")
# remove doublets
kupari003_seurat<-subset(kupari003_seurat, subset= scDblFinder.class == "singlet")

VlnPlot(kupari003_seurat, features = c("nFeature_originalexp", "nCount_originalexp", "subsets_Mito_percent"), ncol = 3)
saveRDS(kupari003_seurat, "../../mouse-nuc-seq/nodose/data/processed/kupari003-seurat-processed-220711.RDS")

#kupari 005 dataset
kupari005<-read10xCounts("kupari2019/10X_18_005/filtered_feature_bc_matrix", col.names = T)
mt_genes<-rownames(kupari005[grep(pattern = "mt-", x = kupari005@rowRanges@elementMetadata$Symbol),])
mt_genes<-mt_genes[-1]
kupari005<-addPerCellQC(kupari005, subset = list(Mito = mt_genes))
plotColData(kupari005, x = "sum", y = "subsets_Mito_percent")
#plotting QC data
x<-plotColData(kupari005, y = "sum")
y<-plotColData(kupari005, y = "detected")
z<-plotColData(kupari005, y = "subsets_Mito_percent")
x + y + z
summary(kupari005$sum)
summary(kupari005$detected)
summary(kupari005$subsets_Mito_percent)
#removing low quality cells 
keep.total<-kupari005$sum >600
keep.top<-kupari005$sum < 75000
keep.n<-kupari005$detected>400
keepmito<-kupari005$subsets_Mito_percent<10
kupari005<-kupari005[,keep.total & keep.n & keepmito & keep.top]
set.seed(1)
kupari005<-scDblFinder(kupari005)
#log normalise to allow for conversion to seurat object 
kupari005<-logNormCounts(kupari005)
kupari005_seurat<-as.Seurat(kupari005)
#add meta data
kupari005_seurat<-AddMetaData(kupari005_seurat, "Kupari", col.name = "Dataset")
kupari005_seurat<-AddMetaData(kupari005_seurat, "Kupari_005", col.name = "Experiment")
kupari005_seurat<-AddMetaData(kupari005_seurat, "5 weeks", col.name = "Age")
kupari005_seurat<-AddMetaData(kupari005_seurat, "Vglut2Cre-Tomato", col.name = "Strain")
# remove doublets
kupari005_seurat<-subset(kupari005_seurat, subset= scDblFinder.class == "singlet")
VlnPlot(kupari005_seurat, features = c("nFeature_originalexp", "nCount_originalexp", "subsets_Mito_percent"), ncol = 3)
saveRDS(kupari005_seurat, "../../mouse-nuc-seq/nodose/data/processed/kupari005-seurat-processed-220711.RDS")


#in house data 
#read in the data, and table with cmo information
inhouse<-read10xCounts("inhouse/filtered_feature_bc_matrix", col.names = T)
cmo<-read.table("inhouse/cell_cmo_assignment_decimal.txt", sep = '\t', header = T, row.names = 3)
#change all NA fields in table to 'unassigned'
cmo[is.na(cmo)] = "unassigned"
#adding meta data to dataset about nutritional condition and left/right nodose
cells.use<-as.character(rownames(cmo))
inhouse<-inhouse[,cells.use]
cmo<-cmo[colnames(inhouse),]
cmo1<-cmo
cmo1<-cbind(cmo1, Nutr.cond = "unassigned")
cmo1<-cmo1[colnames(inhouse),]
index<-0
for (i in cmo1$Sample.Group)
{
  index <- index + 1
  if (i == 1 || i==2) 
  {
    cmo1$Nutr.cond[index] <- "adlib"
  }
  if(i == 3 || i == 4)
  {
    cmo1$Nutr.cond[index] <- "fast"
  }
}
cell.use<-as.character(colnames(inhouse))
cmo1<-cmo1[cell.use,]
inhouse@metadata$Nutr.cond<-NA
inhouse@metadata$Nutr.cond<-cmo1$Nutr.cond

cmo1<-cbind(cmo1, Position = 'unassigned')
cmo1<-cmo1[colnames(inhouse),]
index<-0
for (i in cmo1$Sample.Group)
{
  index <- index + 1
  if (i == 1 || i==3) 
  {
    cmo1$Position[index] <- "left"
  }
  if(i == 2 || i == 4)
  {
    cmo1$Position[index] <- "right"
  }
}
cell.use<-as.character(colnames(inhouse))
cmo1<-cmo1[cell.use,]
inhouse@metadata$Position<-cmo1$Position

#adding mitochondrial QC information
mt_genes<-rownames(inhouse[grep(pattern = "mt-", x = inhouse@rowRanges@elementMetadata$Symbol),])
mt_genes<-mt_genes[-1]
inhouse<-addPerCellQC(inhouse, subset = list(Mito = mt_genes))
plotColData(inhouse, x = "sum", y = "subsets_Mito_percent")
#plotting QC data
x<-plotColData(inhouse, y = "sum")
y<-plotColData(inhouse, y = "detected")
z<-plotColData(inhouse, y = "subsets_Mito_percent")
x + y + z
summary(inhouse$sum)
summary(inhouse$detected)
summary(inhouse$subsets_Mito_percent)
#removing low quality cells 
metakeep<-inhouse@metadata$Nutr.cond!='unassigned'
poskeep<-inhouse@metadata$Position!='unassigned'
keep.total<-inhouse$sum >1000
keep.top<-inhouse$sum < 100000
keep.n<-inhouse$detected>600
keepmito<-inhouse$subsets_Mito_percent<2
inhouse<-inhouse[,keep.total & keep.n & keepmito & keep.top & metakeep & poskeep]
inhouse@metadata$Nutr.cond<-NA
inhouse@metadata$Nutr.cond<-NA

set.seed(1)
inhouse<-scDblFinder(inhouse)
#log normalise to allow for conversion to seurat object 
inhouse<-logNormCounts(inhouse)
inhouse_seurat<-as.Seurat(inhouse)
#add meta data
inhouse_seurat<-AddMetaData(inhouse_seurat, "inhouse", col.name = "Dataset")
inhouse_seurat<-AddMetaData(inhouse_seurat, "inhouse_nucseq", col.name = "Experiment")
inhouse_seurat<-AddMetaData(inhouse_seurat, "7 weeks", col.name = "Age")
inhouse_seurat<-AddMetaData(inhouse_seurat, "WT", col.name = "Strain")
# remove doublets
inhouse_seurat<-subset(inhouse_seurat, subset= scDblFinder.class == "singlet")
VlnPlot(inhouse_seurat, features = c("nFeature_originalexp", "nCount_originalexp", "subsets_Mito_percent"), ncol = 3)

#adding information about left right fed fast
cmo<-read.table("inhouse/cell_cmo_assignment_decimal.txt", sep = '\t', header = T, row.names = 3)
cmo[is.na(cmo)] = "unassigned"
cells.use<-as.character(inhouse_seurat$Barcode)
cmo<-cmo[cells.use,]
cmo<-cmo[rownames(inhouse_seurat@meta.data),]
cmo1<-cmo
cmo1<-cbind(cmo1, Nutr.cond = "unassigned")
cmo1<-cmo1[rownames(inhouse_seurat@meta.data),]
index<-0
for (i in cmo1$Sample.Group)
  {
index <- index + 1
if (i == 1 || i==2) 
  {
cmo1$Nutr.cond[index] <- "adlib"
}
if(i == 3 || i == 4)
  {
cmo1$Nutr.cond[index] <- "fast"
}
}
cell.use<-as.character(inhouse_seurat$Barcode)
cmo1<-cmo1[cell.use,]
identical(rownames(cmo1), inhouse_seurat@assays$originalexp@counts@Dimnames[[2]])
inhouse_seurat<-AddMetaData(inhouse_seurat, cmo1$Nutr.cond, col.name = "Nutr.cond" )

cmo1<-cbind(cmo1, Position = 'unassigned')
cmo1<-cmo1[rownames(inhouse_seurat@meta.data),]
index<-0
for (i in cmo1$Sample.Group)
{
  index <- index + 1
  if (i == 1 || i==3) 
  {
    cmo1$Position[index] <- "left"
  }
  if(i == 2 || i == 4)
  {
    cmo1$Position[index] <- "right"
  }
}
cell.use<-as.character(inhouse_seurat$Barcode)
cmo1<-cmo1[cell.use,]
inhouse_seurat<-AddMetaData(inhouse_seurat, cmo1$Position, col.name = "Position" )
#adding information about groups 
inhouse_seurat<-AddMetaData(inhouse_seurat, cmo$CMO1.., col.name = "pct_cmo1" )
inhouse_seurat<-AddMetaData(inhouse_seurat, cmo$CMO2.., col.name = "pct_cmo2" )
inhouse_seurat<-AddMetaData(inhouse_seurat, cmo$CMO3.., col.name = "pct_cmo3" )
inhouse_seurat<-AddMetaData(inhouse_seurat, cmo$CMO4.., col.name = "pct_cmo4" )
inhouse_seurat<-AddMetaData(inhouse_seurat, cmo$Sample.Group, col.name = "group" )
dim(inhouse_seurat)
saveRDS(inhouse_seurat, "../../mouse-nuc-seq/nodose/data/processed/inhouse-assigned-processed-220712.RDS")

#Buchanan 2022 
buchananleft<-read10xCounts("buchanan2022/left_count/outs/filtered_feature_bc_matrix/", col.names = T)
mt_genes<-rownames(buchananleft[grep(pattern = "mt-", x = buchananleft@rowRanges@elementMetadata$Symbol),])
mt_genes<-mt_genes[-1]
buchananleft<-addPerCellQC(buchananleft, subset = list(Mito = mt_genes))
plotColData(buchananleft, x = "sum", y = "subsets_Mito_percent")
x<-plotColData(buchananleft, y = "sum")
y<-plotColData(buchananleft, y = "detected")
z<-plotColData(buchananleft, y = "subsets_Mito_percent")
x + y + z
summary(buchananleft$sum)
summary(buchananleft$detected)
summary(buchananleft$subsets_Mito_percent)
keep.total<-buchananleft$sum >600
keepsum<-buchananleft$sum<30000
keep.n<-buchananleft$detected>200
keepmito<-buchananleft$subsets_Mito_percent<10
buchananleft<-buchananleft[,keep.total & keep.n & keepmito & keepsum]
set.seed(1)
buchananleft<-scDblFinder(buchananleft)
#log normalise to allow for conversion to seurat object 
buchananleft<-logNormCounts(buchananleft)
buchananleft_seurat<-as.Seurat(buchananleft)
#add meta data
buchananleft_seurat<-AddMetaData(buchananleft_seurat, "Buchanan", col.name = "Dataset")
buchananleft_seurat<-AddMetaData(buchananleft_seurat, "buchanan_left", col.name = "Experiment")
buchananleft_seurat<-AddMetaData(buchananleft_seurat, "6-20 weeks", col.name = "Age")
buchananleft_seurat<-AddMetaData(buchananleft_seurat, "WT", col.name = "Strain")
buchananleft_seurat<-AddMetaData(buchananleft_seurat, "left", col.name = "Position")
# remove doublets
buchananleft_seurat<-subset(buchananleft_seurat, subset= scDblFinder.class == "singlet")
VlnPlot(buchananleft_seurat, features = c("nFeature_originalexp", "nCount_originalexp", "subsets_Mito_percent"), ncol = 3)
saveRDS(buchananleft_seurat, "../../mouse-nuc-seq/nodose/data/processed/buchananleft-processed-220712.RDS")


buchananright<-read10xCounts("buchanan2022/right_count/outs/filtered_feature_bc_matrix/", col.names = T)
mt_genes<-rownames(buchananright[grep(pattern = "mt-", x = buchananright@rowRanges@elementMetadata$Symbol),])
mt_genes<-mt_genes[-1]
buchananright<-addPerCellQC(buchananright, subset = list(Mito = mt_genes))
plotColData(buchananright, x = "sum", y = "subsets_Mito_percent")
x<-plotColData(buchananright, y = "sum")
y<-plotColData(buchananright, y = "detected")
z<-plotColData(buchananright, y = "subsets_Mito_percent")
x + y + z
summary(buchananright$sum)
summary(buchananright$detected)
summary(buchananright$subsets_Mito_percent)
keep.total<-buchananright$sum >600
keepsum<-buchananright$sum<30000
keep.n<-buchananright$detected>200
keepmito<-buchananright$subsets_Mito_percent<10
buchananright<-buchananright[,keep.total & keep.n & keepmito & keepsum]
set.seed(1)
buchananright<-scDblFinder(buchananright)
#log normalise to allow for conversion to seurat object 
buchananright<-logNormCounts(buchananright)
buchananright_seurat<-as.Seurat(buchananright)
#add meta data
buchananright_seurat<-AddMetaData(buchananright_seurat, "Buchanan", col.name = "Dataset")
buchananright_seurat<-AddMetaData(buchananright_seurat, "buchanan_right", col.name = "Experiment")
buchananright_seurat<-AddMetaData(buchananright_seurat, "6-20 weeks", col.name = "Age")
buchananright_seurat<-AddMetaData(buchananright_seurat, "WT", col.name = "Strain")
buchananright_seurat<-AddMetaData(buchananright_seurat, "right", col.name = "Position")
# remove doublets
buchananright_seurat<-subset(buchananright_seurat, subset= scDblFinder.class == "singlet")
VlnPlot(buchananright_seurat, features = c("nFeature_originalexp", "nCount_originalexp", "subsets_Mito_percent"), ncol = 3)
saveRDS(buchananright_seurat, "../../mouse-nuc-seq/nodose/data/processed/buchananright-processed-230401.RDS")

#Zhao 
setwd("zhao2022/")
data<-list()
files<-dir(pattern = "SAMN2460720")
for (i in files)
{
  data[[i]]<-read10xCounts(paste0(print(i),"/outs/filtered_feature_bc_matrix/"), col.names = T)
}
mt_genes<-rownames(data$SAMN24607204[grep(pattern = "mt-", x = data$SAMN24607204@rowRanges@elementMetadata$Symbol),])
mt_genes<-mt_genes[-1]
for (i in files)
{
  data[[i]]<-addPerCellQC(data[[i]], subset = list(Mito = mt_genes))

}
for (i in files)
{
  print(i)
  print("sum")
  print(summary(data[[i]]$sum))
  print("detected")
  print(summary(data[[i]]$detected))
  print("mito")
  print(summary(data[[i]]$subsets_Mito_percent))
}

for (i in files)
{
  keepgenesbottom<-data[[i]]$detected >500
  keepgenestop<-data[[i]]$detected<8000
  keeptranscripts<-data[[i]]$sum>800
  keepmito<-data[[i]]$subsets_Mito_percent<10
  data[[i]]<-data[[i]][,keepgenesbottom & keepgenestop & keepmito & keeptranscripts]
  
}
for (i in files)
{
  set.seed(1)
  data[[i]]<-scDblFinder(data[[i]])
}
#log normalise to allow for conversion to seurat object
data_seurat<-list()
for (i in files)
{
  data[[i]]<-logNormCounts(data[[i]])
  data_seurat[[i]]<-as.Seurat(data[[i]])
}
#add meta data
for (i in files)
{
  data_seurat[[i]]<-AddMetaData(data_seurat[[i]], "Zhao", col.name = "Dataset")
  data_seurat[[i]]<-AddMetaData(data_seurat[[i]], paste0("Zhao_",print(i)), col.name = "Experiment")
  data_seurat[[i]]<-AddMetaData(data_seurat[[i]], "8+ weeks", col.name = "Age")
  data_seurat[[i]]<-AddMetaData(data_seurat[[i]], "WT", col.name = "Strain")
}
# remove doublets
for (i in files)
{
  data_seurat[[i]]<-subset(data_seurat[[i]], subset= scDblFinder.class == "singlet")
}
for (i in files)
  {
  saveRDS(data_seurat[[i]], paste0("../../../mouse-nuc-seq/nodose/data/processed/processed-zhao_", print(i),"220712.RDS"))
}
