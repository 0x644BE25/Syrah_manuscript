######################################################
# MULTI-SPECIES COMPARISON
#
# This assumes that you have downloaded and
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

source('helper_functions.R')

# ================= PARAMS ===========================

datasets <- c('mouse','planarian','chick','human')
versions <- c('original','syrah')
nPCs <- c('mouse'=5,'planarian'=15,'chick'=15,'human'=10)

minUMI <- 50
minUMIs <- c(1,5,10,25,50,100)

ms_colors <- c('mouse'='#6620bb','planarian'='#496efd','chick'='#bb59bf','human'='#ff585f')
ms_colors_lt <- c('mouse'='#b390dd','planarian'='#a4b7fe','chick'='#ddacdf','human'='#ffacaf')

set.seed(61687)

# ================= COLLATE STATS ====================

summary <- NULL

for (ds in datasets) {
  coords <- read.delim(paste0(ds,'_coordinates.tsv'),header=FALSE,row.names=1)
  
  for (vs in versions) {
    if (ds=='planarian') { 
      counts <- getCountsFromNpz(paste0(ds,'_1run_',vs,'_sparse_counts.npz'))
    } else {
      counts <- getCountsFromNpz(paste0(ds,'_',vs,'_sparse_counts.npz'))
    }
    counts <- counts[,colnames(counts) %in% rownames(coords)]
    counts_umis <- colSums(counts)
    
    # compile stats over different min UMI cutoffs
    curr <- do.call(rbind,lapply(minUMIs,\(minUMI){
      beads <- names(counts_umis)[counts_umis>=minUMI]
      if (length(beads)>1) {
        nBeads <- length(beads)
        nUMIs <- sum(counts_umis[beads])
        nGenes <- sum(rowSums(counts[,beads])>0)
        return(data.frame(dataset=ds,version=vs,minUMI,nBeads,nUMIs,nGenes))
      } else {
        nBeads <- 0
        nUMIs <- 0
        nGenes <- 0
        return(data.frame(dataset=ds,version=vs,minUMI,nBeads,nUMIs,nGenes))
      }
    }))
    summary <- rbind(summary,curr)
  }
}
write.csv(summary,'multi-species_comparison_results.csv',row.names=FALSE)

# ================= POISSON TESTS ====================

summary <- read.csv('multi-species_comparison_results.csv')

metrics <- c('nBeads','nUMIs','nGenes')
for (metric in metrics) {
  summary[,paste0('pct_change_',metric)] <- NA
  summary[,paste0('pval_',metric)] <- NA
  summary[,paste0('ci95_low_',metric)] <- NA
  summary[,paste0('ci95_high_',metric)] <- NA
}
for (minUMI in minUMIs) {
  for (ds in datasets) {
    original <- summary[summary$minUMI==minUMI & summary$dataset==ds & summary$version=='original',]
    syrah <- summary[summary$minUMI==minUMI & summary$dataset==ds & summary$version=='syrah',]
    for (metric in c('nBeads','nUMIs','nGenes')) {
      pt <- poisson.test(x=c(syrah[,metric],original[,metric]))
      summary[summary$minUMI==minUMI & summary$dataset==ds & summary$version=='syrah',paste0('pct_change_',metric)] <- 100*(pt$estimate-1)
      summary[summary$minUMI==minUMI & summary$dataset==ds & summary$version=='syrah',paste0('pval_',metric)] <- min_p(pt$p.value)
      summary[summary$minUMI==minUMI & summary$dataset==ds & summary$version=='syrah',paste0('ci95_low_',metric)] <- 100*(pt$conf.int[1]-1)
      summary[summary$minUMI==minUMI & summary$dataset==ds & summary$version=='syrah',paste0('ci95_high_',metric)] <- 100*(pt$conf.int[2]-1)
    }
  }
}
write.csv(summary,'multi-species_comparison_results.csv',row.names=FALSE)

# ================= CREATE SEURAT OBJECTS ============

for (ds in datasets) {
  coords <- read.delim(paste0(ds,'_coordinates.tsv'),header=FALSE,row.names=1)
  pcs_df <- NULL
  
  for (vs in versions) {
    if (ds=='planarian') { 
      counts <- getCountsFromNpz(paste0(ds,'_1run_',vs,'_sparse_counts.npz'))
    } else {
      counts <- getCountsFromNpz(paste0(ds,'_',vs,'_sparse_counts.npz'))
    }
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
