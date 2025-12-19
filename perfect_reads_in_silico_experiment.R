# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289299

######################################################
# PERFECT READS
#
# NOTES: This assumes that you have downloaded the
# 1 run planarian data from GEO accession GSE289299
# (SRA SRX27645703) as well as the bead barcodes and 
# coordinates file from GSE289299 named
# GSE289299_planarian_bead_coordinates.txt.gz
######################################################

# ================= IMPORTS ==========================

# ================= HELPER METHODS ===================

moveIfOld <- function(filename) {
  if (file.exists(filename)) { file.rename(paste0(filename,'.OLD')) }
  return(filename)
}

# ================= PARAMS ===========================

# how many replicates?
nReps <- 5
nReads <- 10^6

# make this the path to the bead barcodes/coordinates file
bcFile <- 'yourPathToCoords.txt.gz'

# starting FASTQ files
r1startFile <- 'yourPathToRead1FASTQ.fastq'
r2startFile<- 'yourPathToRead2FASTQ.fastq'

# YOU SHOULDN'T NEED TO MESS WITH THE FOLLOWING

# computational parameters
batchSize <- 10^5

# linker matching parameters
pattern <- 'JJJJJJJJTCTTCAGCGTTCCCGAGAJJJJJJJTCNNNNNNNNT'
bc1start <- 1
bc1end <- 8
linkstart <- 9
linkend <- 26
bc2start <- 27
bc2end <- 32
TCstart <- 34
TCend <- 35
umistart <- 37
umiend <- 42
linker <- substr(pattern,9,26)
nts <- c('A','C','G','T')

# filtering output FASTQs
r1perfectAll <- moveIfOld('./data/perfect_reads_ALL_R1.fastq')
r2perfectAll <- moveIfOld('./data/perfect_reads_ALL_R2.fastq')

# deduplication output FASTQs
r1perfectDedup <- moveIfOld('./data/perfect_reads_dedup_R1.fastq')
r2perfectDedup <- moveIfOld('./data/perfect_reads_dedup_R2.fastq')

# subsample FASTQs
r1endFiles <- setNames(sapply(1:nReps,function(x){ moveIfOld(paste0('./data/perfect_reads_rep',x,'_perfect','_R1.fastq')) }),1:nReps)
r2endFiles <- setNames(sapply(1:nReps,function(x){ moveIfOld(paste0('./data/perfect_reads_rep',x,'_perfect','_R2.fastq')) }),1:nReps)

# ================= GET PERFECT READS ================
# NOTE: This part exists independently of replicates
# and so is only performed once.

# open connections
conR1 <- file(r1startFile,'r')
conR2 <- file(r2startFile,'r')

# loop
totalReads <- 0
cat('Starting read filtering...\n')
i <- 0
repeat({
  i <- i+1
  r1s <- readLines(conR1,n=batchSize*4)
  r2s <- readLines(conR2,n=batchSize*4)
  if (length(r1s)==0) { break() }
  
  # identify perfect R1s
  idseqs <- sapply(strsplit(r1s[seq(from=1,to=length(r1s),by=4)],' '),\(x){x[2]})
  r1seqs <- r1s[seq(from=2,to=length(r1s),by=4)]
  
  link <- substr(r1seqs,linkstart,linkend)
  bc <- paste0(substr(r1seqs,bc1start,bc1end),substr(r1seqs,bc2start,bc2end))
  tc <- substr(r1seqs,TCstart,TCend)
  lastt <- substr(r1seqs,44,44)
  good_ind <- which(idseqs==demux_to_use & link==linker & tc=='TC' & lastt=='T' & bc %in% barcodes)
  good_ind <- 1+(4*(good_ind-1))
  
  # write good read pairs
  if (length(good_ind)>0) {
    good_inds <- rep(good_ind,each=4)+rep(0:3,length(good_ind))
    whatr1 <- r1s[good_inds]
    whatr2 <- r2s[good_inds]
    if (length(whatr1)!=length(whatr2)) { break() }
    writeR1 <- file(r1perfectAll,'a')
    writeLines(whatr1,writeR1)
    close(writeR1)
    
    writeR2 <- file(r2perfectAll,'a')
    writeLines(whatr2,writeR2)
    close(writeR2)
  }
  # track progess
  totalReads <- totalReads+length(r1seqs)
  if (totalReads%%(10^6)==0) { cat('  ',format(as.POSIXlt(Sys.time())),' ',totalReads,'reads processed\n') }
  
  # end condition
  if (length(r1seqs)<batchSize) { break() }
})
closeAllConnections()

# ================= DEDUPLICATE PERFECT READS ========
# NOTE: This part ALSO exists independently of 
# replicates and so is only performed once.

# open file connections
conR1 <- file(r1perfectAll,'r')
conR2 <- file(r2perfectAll,'r')

# loop
totalReads <- 0
cat('Starting read deduplication...\n')
i <- 0
seen_bcumis <- NULL
repeat({
  i <- i+1
  r1s <- readLines(conR1,n=batchSize*4)
  r2s <- readLines(conR2,n=batchSize*4)
  if (length(r1s)==0) { break() }
  
  # identify unique R1s
  r1seqs <- r1s[seq(from=2,to=length(r1s),by=4)]
  
  bcumis <- paste0(substr(r1seqs,bc1start,bc1end),substr(r1seqs,bc2start,bc2end),substr(r1seqs,umistart,umiend))
  good_ind <- match(setdiff(bcumis,seen_bcumis), bcumis)
  seen_bcumis <- c(seen_bcumis,bcumis[good_ind])
  good_ind <- 1+(4*(good_ind-1))
  
  # write unique read pairs
  if (length(good_ind)>0) {
    good_inds <- rep(good_ind,each=4)+rep(0:3,length(good_ind))
    whatr1 <- r1s[good_inds]
    whatr2 <- r2s[good_inds]
    if (length(whatr1)!=length(whatr2)) { break() }
    writeR1 <- file(r1perfectDedup,'a')
    writeLines(whatr1,writeR1)
    close(writeR1)
    
    writeR2 <- file(r2perfectDedup,'a')
    writeLines(whatr2,writeR2)
    close(writeR2)
  }
  
  # track progress
  totalReads <- totalReads+length(r1seqs)
  if (totalReads%%(10^6)==0) { cat('  ',format(as.POSIXlt(Sys.time())),' ',totalReads,'reads processed\n') }
  
  # end condition(s)
  if (length(r1seqs)<batchSize) { break() }
})
closeAllConnections()

# ================= SUBSAMPLE ========================
# deduped read FASTQs are small enough that we can
# read in the whole thing and iterate over subsampling

r1 <- readLines(r1perfectDedup)
r2 <- readLines(r2perfectDedup)

inds <- seq(from=1,to=length(r1),by=4)
chosen_inds <- replicate(nReps,sort(sample(inds,nReads,replace=FALSE)))


r1 <- readLines('./data/fastqs/H5_perfect_unique_r1.fastq')
r2 <- readLines('./data/fastqs/H5_perfect_unique_r2.fastq')

inds <- seq(from=1,to=length(r1),by=4)
chosen_inds <- replicate(length(reps),sort(sample(inds,nreads,replace=FALSE)))

for (n in reps) {
  cat(n)
  rep_name <- reps[n]
  cat(rep_name)
  rep_inds <- chosen_inds[,n]
  rep_lines <- rep(rep_inds,each=4)+rep(0:3,length(rep_inds))
  writeLines(r1[rep_lines],paste0('./data/fastqs/H51M-',rep_name,'_perfect_r1.fastq'))
  writeLines(r2[rep_lines],paste0('./data/fastqs/H51M-',rep_name,'_perfect_r2.fastq'))
}

# ================= SUBSET +ADD DELETIONS ====================

for (n in reps) {
  
  read1fastq <- paste0('./data/fastqs/H51M-',n,'_perfect_r1.fastq')
  r1outs <- setNames(paste0(write_dir,'H51M-',n,'_',c(paste0('del',delrates),'nonsense'),'_r1.fastq'),c(delrates,'nonsense'))
  for (x in r1outs) { if (file.exists(x)) { file.rename(x,paste0(x,'.OLD')) }}
  
  len <- nchar(readLines(read1fastq,n=2)[2])
  ts <- rep('T',len)
  
  r1s <- readLines(read1fastq)
  seq_ind <- seq(from=2,to=length(r1s),by=4)
  seqs <- r1s[seq_ind]
  
  for (dr in delrates) {
    dels <- replicate(length(seqs),(sample(c(TRUE,FALSE),size=len,replace=TRUE,prob=c(1-dr,dr))))
    ind <- which(apply(dels,2,mean)<1)
    w_del <- sapply(ind,\(i){ paste0(c(strsplit(seqs[i],'')[[1]][dels[,i]],ts)[1:len],collapse='') })
    seqs[ind] <- w_del
    curr <- r1s
    curr[seq_ind] <- seqs
    writeR1 <- file(r1outs[as.character(dr)],'a')
    writeLines(curr,writeR1)
    close(writeR1)
  }
  dr <- 'nonsense'
  rands <- replicate(length(seqs),paste0(sample(nts,size=len,replace=TRUE),collapse=''))
  curr[seq_ind] <- rands
  writeLines(curr,r1outs[as.character(dr)])
  closeAllConnections()
}

# After this you need to run Syrah on each of the read 1 files
# and then run everything through the WARP pipeline using
# the SchmedS3 custom reference available via 
# gcloud storage cp 'gs://fc-74bc109e-8b3b-4876-a002-ba0ffa3fb758/uploads/planarian/modified_star2.7.10a-Schmidtea_mediterranea-GENCODE-build-smed_20140614-schmiddtea_mediterranea.tar'