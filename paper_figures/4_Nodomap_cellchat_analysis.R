# ============================================
# Script: Nodomap_cellchat_analysis.R
# Purpose: Run cellchat analysis on combined nodomap and mouse hindbrain dataset to identify possible communication between nodose ganglion neurons and hindbrain neurons. 
# Inputs:  - mouse-hindbrain-fedfast-2023.RDS (mouse hindbrain dataset)
#          - diet_nodomap.RDS (nodomap dataset)
#          - mouse_genes_names_Ens100_rmdup.txt (gene_to_ens matching table)
# ============================================

# Load required packages ----
library(Seurat)
library(ggplot2)
library(CellChat)

# Set up paths and config ----
data_dir = ("/cephfs/scratch/bruening_scratch/gdowsett/data/nodomap/")

# Load data ----
hb<-readRDS(paste0(data_dir, 'mouse-hindbrain-fedfast-2023.RDS'))
hb <- UpdateSeuratObject(hb)

nodomap<-readRDS(paste0(data_dir, 'diet_nodomap_cleaned.RDS'))
nodomap <- UpdateSeuratObject(nodomap)

gene_names1<-read.delim(paste0(data_dir, "mouse_genes_names_Ens100_rmdup.txt"),row.names=2)
gene_names3<-read.delim(paste0(data_dir, "mouse_genes_names_Ens100_rmdup.txt"),row.names=1)
gene_names<-read.delim(paste0(data_dir, "mouse_genes_names_Ens100_rmdup.txt"),row.names=3)

# internal functions ----
genes_to_ens<-function(gene_symbol){
  return(as.character(gene_names[gene_symbol,1]))
}

calc_percent_cluster <- function(seurat_object, genes, metadata_name, assay = NULL, layer = "data") {
  # Pick default assay if not provided
  if (is.null(assay)) {
    assay <- DefaultAssay(seurat_object)
  }
  
  # Fetch expression values for selected genes (per cell)
  gene_expr <- FetchData(seurat_object, vars = genes, layer = layer)
  
  # Extract cluster info
  clusters <- seurat_object@meta.data[, metadata_name]
  
  # Average expression per cluster (this is *not* pseudobulk)
  avex_percluster <- AverageExpression(
    seurat_object,
    features = genes,
    group.by = metadata_name,
    assays = assay,
    slot = layer
  )
  
  # AverageExpression returns a list → extract correct assay
  avex_percluster <- data.frame(t(avex_percluster[[assay]]))
  
  # Round and rename
  colnames(avex_percluster) <- paste0("AvEx_", colnames(avex_percluster))
  avex_percluster <- round(avex_percluster, 2)
  
  # Convert expression to binary (expressed vs not)
  gene_expr[gene_expr > 0] <- 1
  gene_expr[gene_expr < 1] <- 0
  
  # Count total cells per cluster
  total_cells_per_cluster <- table(clusters)
  total_cells_per_cluster <- subset(total_cells_per_cluster, total_cells_per_cluster > 0)
  
  # Count expressing cells per cluster
  cells_expressing_gene_per_cluster <- rowsum(gene_expr, clusters)
  
  # Percent expressing
  percent_expressing_gene_per_cluster <- cells_expressing_gene_per_cluster / total_cells_per_cluster * 100
  percent_expressing_gene_per_cluster <- round(percent_expressing_gene_per_cluster, 2)
  
  # Combine results
  result <- cbind(percent_expressing_gene_per_cluster, avex_percluster)
  
  return(as.data.frame(result))
}

convert_colnames <- function(df) {
  new_colnames <- sapply(colnames(df), function(ensid) {
    gene_symbol <- gene_names1[gene_names1$Gene.stable.ID == ensid, 'Gene.name']
    if (length(gene_symbol) == 0) {
      return(ensid)  # Return original ENSID if no matching gene symbol
    } else {
      return(as.character(gene_symbol))
    }
  })
  colnames(df) <- new_colnames
  return(df)
}


# Analysis ----

# Step 1: Prepare hindbrain object for CellChat
Idents(hb)<-'cell_type'
hb<-subset(hb, idents = 'Neurons')
hb_norm<-GetAssayData(hb, assay = 'RNA', layer = 'data')
#print first  IDs for downstream sanity check 
print(rownames(hb_norm)[1:20])
#Change ENSID to gene names in the matrix
match_indices <- match(rownames(hb_norm), gene_names1$Gene.stable.ID)
genes <- ifelse(!is.na(match_indices), gene_names1$Gene.name[match_indices], rownames(hb_norm))
rownames(hb_norm)<-genes
#print new gene names for sanity check 
print(rownames(hb_norm)[1:20])
#check some names are correct
genes_to_ens('Sox17')
genes_to_ens('Xkr4')
#get cluster names
Idents(hb)<-'cluster_names'
labels<-Idents(hb)
#Make metadata
meta<-data.frame(group = labels, row.names = names(labels))
#prefix cluster names with HB
meta$group <- paste("HB_", meta$group, sep = "")

# Step 2: Pepare Nodose object
Idents(nodomap)<-'cell_identity'
nodomap<-subset(nodomap, idents = 'Nodose ganglion neuron')
nodo_norm<-GetAssayData(nodomap, assay = 'originalexp', slot = 'data')
#print first  IDs for downstream sanity check 
print(rownames(nodo_norm)[1:20])
#Change ENSID to gene names in the matrix
match_indices <- match(rownames(nodo_norm), gene_names1$Gene.stable.ID)
genes <- ifelse(!is.na(match_indices), gene_names1$Gene.name[match_indices], rownames(nodo_norm))
rownames(nodo_norm)<-genes
#print new gene names for sanity check 
print(rownames(nodo_norm)[1:20])
#check some names are correct
genes_to_ens('Sox17')
genes_to_ens('Xkr4')
#get cluster names
Idents(nodomap)<-'cluster_number'
labels<-Idents(nodomap)
#Make metadata
meta_nodo<-data.frame(group = labels, row.names = names(labels))
meta_nodo$group<-paste("Nodo_", meta_nodo$group, sep = "")

# Step 3: Create merged normalised object
common_genes <- intersect(rownames(hb_norm), rownames(nodo_norm))
hb_norm <- hb_norm[common_genes, ]
nodo_norm <- nodo_norm[common_genes, ]
#check gene orders are the same 
identical(rownames(hb_norm), rownames(nodo_norm))
merged <- RowMergeSparseMatrices(hb_norm, nodo_norm)
# Ensure metadata rows match the cell names in the matrices
meta_hb <- meta[colnames(hb_norm), , drop = FALSE]
meta_nodo <- meta_nodo[colnames(nodo_norm), , drop = FALSE]
meta_append <- rbind(meta_hb, meta_nodo)
all(rownames(meta_append) == colnames(merged))

# Step 4: Run cellchat
cellchat<-createCellChat(merged, meta = meta_append, group.by = 'group')
#Use the non-protein signalling pathways only for the database (this is generally neuronal and metabolic signalling)
CellChatDB<-CellChatDB.mouse
non_protein_secreted<-subsetDB(CellChatDB, search = c('Non-protein Signaling', 'Secreted Signaling'), key = 'annotation')
cellchat@DB<-CellChatDB
#subset dataset by only genes in the non_protein signalling
cellchat<-subsetData(cellchat)
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
#Inference of cell-cell communication network 
#using triMean as this produces fewer but stronger interactions 
#triMean says that the average gene expression is 0 if the % of cells in 1 group expressin a gene is less than 25%
cellchat <- computeCommunProb(cellchat, type = "triMean")
#Filter the cell-cell communication (default number of cells is 10)
cellchat <- filterCommunication(cellchat, min.cells = 20)
#infer cell-cell communication at a signaling pathway level 
cellchat<-computeCommunProbPathway(cellchat)
#calculate an aggregated cell-cell communication network
cellchat<-aggregateNet(cellchat)

# Step 5: Save out cell chat object 
saveRDS(cellchat, paste0(data_dir, 'nodomap_hb_cellchat_object.RDS'))
#cellchat <-readRDS(paste0(data_dir, 'nodomap_hb_cellchat_object.RDS'))

# Step 6: Look into output 
clusters <- unique(cellchat@meta$group) 
nodo_clusters <- clusters[startsWith(clusters, 'Nodo_')]
hb_clusters <- clusters[startsWith(clusters, 'HB_')]

cellchat@netP[['pathways']]
extractEnrichedLR(cellchat, signaling = cellchat@netP[['pathways']], geneLR.return = TRUE)
netAnalysis_contribution(cellchat, signaling = c(cellchat@netP[['pathways']]), title = 'Contribution of each LR pairs')
extractEnrichedLR(cellchat, signaling = 'NRXN', geneLR.return = FALSE)

# Step 7: Group pathways depending on frequency
interaction_matrix <- subsetCommunication(cellchat, slot.name = 'netP', sources.use = nodo_clusters, targets.use = hb_clusters)
pathways<-unique(interaction_matrix$pathway_name)
#Initialise empty data frame
interaction_classification<-data.frame(
  Pathway = character(),
  Source_num = numeric(),
  Source_class = character(),
  Target_num = numeric(),
  Target_class = character())


for (i in pathways){
  x<-subset(interaction_matrix, interaction_matrix$pathway_name == i)
  snum<-length(unique(x$source))
  if (snum > length(unique(nodo_clusters))*0.25){
    sclass <- 'many'
  }
  else {
    sclass <- 'few'
  }
  tnum<-length(unique(x$target))
  if (tnum > length(unique(hb_clusters))*0.25){
    tclass <- 'many'
  }
  else {
    tclass <- 'few'
  }
  new_row<-data.frame(
    Pathway = i,
    Source_num = snum,
    Source_class = sclass,
    Target_num = tnum,
    Target_class = tclass)
  interaction_classification<-rbind(interaction_classification, new_row)
}

# Step 8: Rebuttal work: 
# First we want tolook at Oxtr+ nodose clusters to PPG hindbrain neurons 
oxtr <- calc_percent_cluster(nodomap, genes = genes_to_ens('Oxtr'), metadata_name = 'cluster_number')

# top 4 clusters have expr >10% 
top4 <- head(rownames(oxtr)[order(oxtr[, 2], decreasing = TRUE)], 4)
top4 <- paste0("Nodo_", sub("\\..*$", "", top4))
# look at interactions between top 4 and ppg neurons 
oxtr_to_ppg <- interaction_matrix[interaction_matrix$source %in% top4 & interaction_matrix$target == "HB_NE_Gcg/Prlr", ]

oxtr_to_ppg <- oxtr_to_ppg %>%
  left_join(
    interaction_classification %>% select(Pathway, Source_class, Target_class),
    by = c("pathway_name" = "Pathway")
  )
x <- subsetCommunication(cellchat, slot.name = 'net', sources.use = nodo_clusters, targets.use = hb_clusters)
oxtr_to_ppg_annot <- oxtr_to_ppg %>%
  left_join(
    x %>% select(pathway_name, annotation) %>% distinct(),
    by = "pathway_name"
  )
oxtr_to_ppg_sub <- oxtr_to_ppg_annot %>%
  filter(annotation %in% c("Secreted Signaling", "Non-protein Signaling"))

write.csv(oxtr_to_ppg, paste0(data_dir, 'Oxtr_to_ppg_cellchat.csv'), row.names= F)

# Step 9: Plotting for supplementary figure:
# Define sources and targets
tac_source <- c('Nodo_NGN10', 'Nodo_NGN12', 'Nodo_NGN14', 'Nodo_NGN16')
tac_target <- c('HB_NE_Gcg/Prlr', 'HB_NE_Grin2c/Osbpl3', 'NB_NE_Lhx2/Lhx9', 'HB_NE_Mylk/Lpar1')

# pick your groups (sources + targets)
nodes.keep <- unique(c(tac_source, tac_target))

# subset the CellChat object to only those identities
cellchat.sub <- subsetCellChat(cellchat, idents.use = nodes.keep)

pdf(paste0(data_dir,"TAC_circle_subset.pdf"), width = 8, height = 8)
netVisual_aggregate(
  cellchat.sub,
  sources.use = tac_source,
  targets.use = tac_target,
  signaling = "TAC",
  layout = "circle"
)
dev.off()

enho_target <- unique(as.character(interaction_matrix[interaction_matrix$pathway_name == 'ENHO', "target"]))
enho_source <- unique(as.character(interaction_matrix[interaction_matrix$pathway_name == 'ENHO', "source"]))

nodes.keep <- unique(c(enho_source, enho_target))

# subset the CellChat object to only those identities
cellchat.sub <- subsetCellChat(cellchat, idents.use = nodes.keep)

pdf(paste0(data_dir,"ENHO_circle_subset.pdf"), width = 8, height = 8)
netVisual_aggregate(
  cellchat.sub,
  sources.use = enho_source,
  targets.use = enho_target,
  signaling = "ENHO",
  layout = "circle"
)
dev.off()

# Save output ----
x <- subsetCommunication(cellchat, slot.name = 'net', sources.use = nodo_clusters, targets.use = hb_clusters)
write.table(x, paste0(data_dir, 'cellchat-table.txt'), sep = '\t', row.names = F)

