######################################################
# DE ANALYSIS AND CELL TYPE PREDICTIONS
#
# This assumes that you have downloaded and
# unpacked the WARP_output planarian data files into 
# your working folder as well as the 
# helper_functions.R script
#
# All are available from:
# https://github.com/0x644BE25/Syrah_manuscript
#
# The LabelTransfer portion requires a prepared Seurat
# object of the pertinent Benham-Pyle et al. 2021
# metadata-ed, aligned to schMedS3, and SCTransformed.
# This file is way too big for GitHub, but easily
# available by request to cbrewster@stowers.org
#
# If you don't have that data, it will use the 
# precomputed labeltransfer results included in
# the WARP_output planarian.tar.gz data.
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


tissueMarkers <- c('piwi-1'='h1SMcG0013999','mat2a'='h1SMcG0005079')

knownGenes <- list('piwi-1'=c('h1SMcG0002117',
                              'h1SMcG0003097',
                              'h1SMcG0006433',
                              'h1SMcG0007496',
                              'h1SMcG0008035',
                              'h1SMcG0009590',
                              'h1SMcG0010823',
                              'h1SMcG0013627',
                              'h1SMcG0013746',
                              'h1SMcG0013999',
                              'h1SMcG0014251',
                              'h1SMcG0014692',
                              'h1SMcG0015722',
                              'h1SMcG0016303',
                              'h1SMcG0018886',
                              'h1SMcG0019482',
                              'h1SMcG0020811',
                              'h1SMcG0020880',
                              'h1SMcG0021424',
                              'h1SMcG0021692',
                              'h1SMcG0021969',
                              'h1SMcG0022103'),
                   'mat2a'=c('h1SMcG0001385',
                             'h1SMcG0002222',
                             'h1SMcG0002269',
                             'h1SMcG0004299',
                             'h1SMcG0007332',
                             'h1SMcG0008011',
                             'h1SMcG0009175',
                             'h1SMcG0011053',
                             'h1SMcG0011829',
                             'h1SMcG0013627',
                             'h1SMcG0015147',
                             'h1SMcG0015360',
                             'h1SMcG0017763',
                             'h1SMcG0021980'))

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
write.csv(summary,'additional_sequencing_results.csv')

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

# ================= DIFFERENTIAL EXPRESSION ==========

# by cluster
markers <- NULL
for (sq in sequencing) {
  for (vs in versions) {
    seu <- readRDS(paste0('planarian_',sq,'_',vs,'_seurObj.rds'))
    seu@active.ident <- seu$seurat_clusters
    seu@active.assay <- 'SCT'
    m <- FindAllMarkers(seu)
    m <- m[m$p_val_adj<.05,]
    m$gene <- rownames(m)
    m$sequencing <- sq
    m$version <- vs
    markers <- rbind(markers,m)
  }
}
write.csv(markers,'planarian_cluster_marker_results.csv',row.names=FALSE)

# known tissues
markers <- NULL
for (n in names(tissueMarkers)) {  
  gene <- tissueMarkers[n]
  known <- knownGenes[[n]]
  
  for (sq in sequencing) {
    for (vs in versions) {
      seu <- readRDS(paste0('planarian_',sq,'_',vs,'_seurObj.rds'))
      pos_beads <- Cells(seu)[seu[['RNA']]$counts[gene,]>0]
      neg_beads <- Cells(seu)[seu[['RNA']]$counts[gene,]==0]
      m <- FindMarkers(seu,ident.1=pos_beads,ident.2=neg_beads,only.pos=TRUE)    
      m$gene <- rownames(m)
      m <- m[m$gene %in% known & m$p_val_adj<.05,]
      m$sequencing <- sq
      m$version <- vs
      m$marker_name <- n
      m$marker_id <- gene
      markers <- rbind(markers,m)
    }
  }
}
write.csv(markers,'known_tissue_marker_DE_results.csv',row.names=FALSE)

# ================= CELL TYPE PREDICTION =============

# do predictions
if (file.exists('TRACS_UI_seurObj.rds')) {
ui <- readRDS('TRACS_UI.rds') # YOU DON'T HAVE THIS
  
  lt_scores <- NULL
  for (sq in sequencing) {
    for (vs in versions) {
      seu <- readRDS(paste0('planarian_',sq,'_',vs,'_seurObj.rds'))
      anchors <- FindTransferAnchors(reference=ui,query=seu,normalization.method='SCT')
      lt <- TransferData(anchors,ui$Tissue)
      df <- data.frame(sequencing=sq,version=vs,barcode=rownames(lt),prediction=lt$predicted.id,score=lt$prediction.score.max)
      lt_scores <- rbind(lt_scores,df)
    }
  }
  write.csv(lt_scores,'./cell_type_prediction.scores')
} else {
  message('You do not have required data to perform LabelTransfer.\nUsing pre-computed data.\nSee head of this script for details.')
}

# stats
lt_scores <- read.csv('./planarian_labeltransfer_results.csv')
ctrl <- lt_scores[lt_scores$sequencing=='1run' & lt_scores$version=='original','score']
summary <- NULL
for (sq in sequencing) {
  for (vs in versions) {
    curr <- lt_scores[lt_scores$sequencing==sq & lt_scores$version==vs,'score']
    tt <- t.test(curr,ctrl)
    df <- data.frame(sequencing=sq,version=vs,
                      pct_change=unname(pct_change(tt$estimate[1],tt$estimate[2])),
                      pval=min_p(tt$p.value),
                      ci95_low=unname(pct_change(tt$estimate[2]+tt$conf.int[1],tt$estimate[2])),
                      ci95_high=unname(pct_change(tt$estimate[2]+tt$conf.int[2],tt$estimate[2])))
    summary <- rbind(summary,df)
  }
}
write.csv(summary,'planarian_cell_type_prediction_summary.csv',row.names=FALSE)
