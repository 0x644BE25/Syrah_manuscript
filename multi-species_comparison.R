######################################################
# MULTI-SPECIES COMPARISON
#
# GOAL: This assumes that you have downloaded and
# unpacked the WARP_output human, mouse, chick, and
# planarian data files into your working folder as 
# well as the helper_functions.R script. 
#
# All are available from:
# https://github.com/0x644BE25/Syrah_manuscript
######################################################

# ================= IMPORTS ==========================

if (!require('Seurat')) { install.packages('Seurat') }
if (!require('ggplot2')) { install.packages('ggplot2') }

source(helper_functions.R)

# ================= PARAMS ===========================

datasets <- c('mouse','planarian','chick','human')
versions <- c('original','syrah')
nPCs <- c('mouse'=5,'planarian'=15,'chick'=15,'human'=10)

minUMI <- 50

ms_colors <- c('mouse'='#6620bb','planarian'='#496efd','chick'='#bb59bf','human'='#ff585f')
ms_colors_lt <- c('mouse'='#b390dd','planarian'='#a4b7fe','chick'='#ddacdf','human'='#ffacaf')

set.seed(61687)

# ================= ANALYSIS =========================

for (ds in datasets) {
  coords <- read.delim(paste0(ds,'_coordinates.tsv'),header=FALSE,row.names=1)
  pcs_df <- NULL
  
  for (vs in versions) {
    counts <- getCountsFromNpz(paste0(ds,'_',vs,'_sparse_counts.npz'))
    counts <- counts[,colnames(counts) %in% rownames(coords)]
    counts <- counts[,colSums(counts)>=minUMI]
    
    # standard single-cell processing
    seu <- CreateSeuratObject(counts=counts,project=paste0(ds,'_',vs))
    seu <- SCTransform(seu)
    seu <- RunPCA(seu,npcs=nPCs[ds])
    seu <- RunUMAP(seu,dims=1:nPCs[ds],seed.use=61687)
    seu <- FindNeighbors(seu,dims=1:nPCs[ds])
    seu <- FindClusters(seu)
    
    # add spatial embedding
    puck <- as.matrix(coords[Cells(seu),]); colnames(coords) <- c('PUCK_1','PUCK_2')
    seu[['puck']] <- CreateDimReducObject(embeddings=puck,key='PUCK_',assay='RNA')
    
    # save Seurat object
    saveRDS(seu,paste0(ds,'_',vs,'_seurObj.rds'),compress=FALSE)
    
    # collate PC info
    pcs_df <- rbind(pcs_df,data.frame(version=vs,PC=1:nPCs[ds],stdev=seu[['pca']]@stdev))
  }
  write.csv(pcs_df,paste0(ds,'_PC_stdevs.csv'))
}
