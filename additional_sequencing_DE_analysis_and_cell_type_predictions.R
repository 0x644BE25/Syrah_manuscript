######################################################
# DE ANALYSIS AND CELL TYPE PREDICTIONS
#
# This assumes that you have downloaded and
# unpacked the WARP_output planarian data files into 
# your working folder as well as the 
# helper_functions.R script. 
#
# All are available from:
# https://github.com/0x644BE25/Syrah_manuscript
######################################################

# ================= IMPORTS ==========================

if (!require('Seurat')) { install.packages('Seurat') }
if (!require('ggplot2')) { install.packages('ggplot2') }

source('helper_functions.R')

# ================= PARAMS ===========================

sequencing <- c('1run','2runs')
versions <- c('original','syrah')
nPCs <- 15

minUMI <- 50
minUMIs <- c(1,5,10,25,50,100)

set.seed(61687)

coords <- read.delim('planarian_coordinates.tsv',header=FALSE,row.names=1)

# ================= COLLATE STATS ====================

summary <- NULL
for (sq in sequencing) {
  for (vs in versions) {
    counts <- getCountsFromNpz(paste0('planarian_',sq,'_',vs,'_sparse_counts.npz'))
    counts <- counts[,colnames(counts) %in% rownames(coords)]
    counts_umis <- colSums(counts)
    
    # compile stats over different min UMI cutoffs
    curr <- do.call(rbind,lapply(minUMIs,\(minUMI){
      beads <- names(counts_umis)[counts_umis>=minUMI]
      if (length(beads)>1) {
        nBeads <- length(beads)
        nUMIs <- sum(counts_umis[beads])
        nGenes <- sum(rowSums(counts[,beads])>0)
        return(data.frame(sequencing=sq,version=vs,minUMI,nBeads,nUMIs,nGenes))
      } else {
        nBeads <- 0
        nUMIs <- 0
        nGenes <- 0
        return(data.frame(sequencing=sq,version=vs,minUMI,nBeads,nUMIs,nGenes))
      }
    }))
    summary <- rbind(summary,curr)
  }
}
write.csv(summary,'additional_sequencing_results.csv',row.names=FALSE)

# ================= POISSON TESTS ====================

summary <- read.csv('additional_sequencing_results.csv')

metrics <- c('nBeads','nUMIs','nGenes')
for (metric in metrics) {
  summary[,paste0('pct_change_',metric)] <- 0
  summary[,paste0('pval_',metric)] <- 0
  summary[,paste0('ci95_low_',metric)] <- 0
  summary[,paste0('ci95_high_',metric)] <- 0
}
for (minUMI in minUMIs) {
  control <- summary[summary$sequencing=='1run' & summary$version=='original' & summary$minUMI==minUMI,]
  for (sq in sequencing) {
    for (vs in versions) {
      test <- summary[summary$sequencing==sq & summary$version==vs & summary$minUMI==minUMI,]
      for (metric in c('nBeads','nUMIs','nGenes')) {
        pt <- poisson.test(x=c(test[,metric]),control[,metric])
        summary[summary$sequencing==sq & summary$version==vs & summary$minUMI==minUMI,paste0('pct_change_',metric)] <- 100*(pt$estimate-1)
        summary[summary$sequencing==sq & summary$version==vs & summary$minUMI==minUMI,paste0('pval_',metric)] <- pt$p.value
        summary[summary$sequencing==sq & summary$version==vs & summary$minUMI==minUMI,paste0('ci95_low_',metric)] <- 100*(pt$conf.int[1]-1)
        summary[summary$sequencing==sq & summary$version==vs & summary$minUMI==minUMI,paste0('ci95_high_',metric)] <- 100*(pt$conf.int[2]-1)
      }
    }
  }
}
write.csv(summary,'./additional_sequencing_results.csv')

# ================= CREATE SEURAT OBJECTS ============

coords <- read.delim('planarian_coordinates.tsv',header=FALSE,row.names=1)
colnames(coords) <- c('PUCK_1','PUCK_2')

for (sq in sequencing) {
  for (vs in versions) {
    counts <- getCountsFromNpz(paste0('planarian_',sq,'_',vs,'_sparse_counts.npz'))
    counts <- counts[,colnames(counts) %in% rownames(coords)]
    counts <- counts[,colSums(counts)>=minUMI]
    
    # standard single-cell processing
    seu <- CreateSeuratObject(counts=counts,project=paste0('planarian_',sq,'_',vs))
    seu <- SCTransform(seu)
    seu <- RunPCA(seu,npcs=nPCs)
    seu <- RunUMAP(seu,dims=1:nPCs,seed.use=61687)
    seu <- FindNeighbors(seu,dims=1:nPCs)
    seu <- FindClusters(seu)
    
    # add spatial embedding
    puck <- as.matrix(coords[Cells(seu),]); colnames(coords) <- c('PUCK_1','PUCK_2')
    seu[['puck']] <- CreateDimReducObject(embeddings=puck,key='PUCK_',assay='RNA')
    
    # save Seurat object
    saveRDS(seu,paste0('planarian_',sq,'_',vs,'_seurObj.rds'),compress=FALSE)
  }
}
