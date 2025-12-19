######################################################
# FIRST THOUGHTS
#
# GOAL: We want to generate an artificial "slide-seq"
# dataset from an image. 
#   * Use image agt appropriate resolution based
#     on roughly 30000 beads in a standard slide-seq
#     dataset = diameter of ~300 beads = 300x300 px
#  * Use pixel locations for bead coordinates with 
#    randomized barcodes. 
#  * Use RGB values for gene expression.
#  * Probabalistically generate reads, do deletions
#    and backfill with 'T', then randomly subsample 
#    to ~10%
######################################################

# ================= IMPORTS ==========================

library(doParallel)
library(foreach)
library(ggplot2)
library(png)

# ================= PARAMS ===========================

images <- c('fruit4','frenchflag','pearlearring','monalisa')
nts <- c('A','C','G','T')
puckdiameter <- 5500

# ================= METHODS ==========================

hexColor <- function(r,g,b) {
  r <- stringr::str_pad(as.character(as.hexmode(ceiling(1*r))),width=2,side='left',pad='0')
  g <- stringr::str_pad(as.character(as.hexmode(ceiling(1*g))),width=2,side='left',pad='0')
  b <- stringr::str_pad(as.character(as.hexmode(ceiling(1*b))),width=2,side='left',pad='0')
  return(paste0('#',r,g,b))
}

# ================= INIT DATA ========================

# ================= ANALYSIS =========================

for (imgName in images) {
  
  # process images ===============
  img <- readPNG(paste0('./data/',imgName,'_300.png'))
  coords <- data.frame(y=rep(ncol(img):1,nrow(img)),x=rep(1:nrow(img),each=ncol(img)))
  bcs <- setdiff(unique(sapply(1:(2*nrow(coords)),\(x){ paste0(sample(nts,14,replace=TRUE),collapse='') })),'TTTTTTTTTTTTTT')
  rownames(coords) <- bcs[1:nrow(coords)]
  rownames(coords)[coords$x==275 & coords$y==275] <- 'TTTTTTTTTTTTTT'
  coords$R <- 255*c(img[,,1])
  coords$G <- 255*c(img[,,2])
  coords$B <- 255*c(img[,,3])
  coords$hex <- hexColor(coords$R,coords$G,coords$B)
  coords$nUMI <- coords$R+coords$G+coords$B
  p1 <- ggplot(coords,aes(x=x,y=y,color=I(hex))) + 
    geom_point(size=1) + 
    coord_fixed() +
    theme_void()
  coeff <- puckdiameter/max(coords$x)
  coords$puck_x <- coords$x*coeff
  coords$puck_y <- coords$y*coeff
  coords <- coords[sample(1:nrow(coords)),]
  write.csv(coords,paste0('./data/',imgName,'_coords_colors_bcs.csv'),row.names=TRUE)
  write.table(coords[,c('puck_x','puck_y')],paste0('./data/',imgName,'_coords.tsv'),sep='\t',row.names=TRUE,col.names=FALSE,quote=FALSE)
  
  # generate reads =================
  cl <- makeForkCluster(5)
  registerDoParallel(cl)
  x <- foreach(k=1:45) %dopar% {
    nts <- c('A','C','G','T')
    baseid <- paste0('@H3GLCDSX3:1:1101:',10040+k,':')
    linker <- 'TCTTCAGCGTTCCCGAGA'
    dr <- .01
    cols <- read.csv('./data/mmus_seqs.csv',header=FALSE)
    cols <- c('R'=cols[cols$V1=='R','V3'][1],
              'G'=cols[cols$V1=='G','V3'][1],
              'B'=cols[cols$V1=='B','V3'][1])
    
    r1qual <- paste0(rep('D',42),collapse='')
    r2qual <- paste0(rep('D',nchar(cols['R'])),collapse='')
    ts <- rep('T',42)
    coords <- read.csv(paste0('./data/',imgName,'_coords_colors_bcs.csv'),row.names=1)
    cat(k,nrow(coords),'\n')
    start <- 2000*(k-1)+1
    end <- 2000*k
    coords <- coords[start:end,]
    coords <- coords[coords$nUMI>0,]
    
    cat(start,end,'\n')
    for (i in 1:nrow(coords)) {
      currow <- coords[i,]
      bc <- rownames(currow)
      nUMI <- currow$nUMI
      if (nUMI>0) {
        r2s <- c(rep(cols['R'],currow['R']),
                 rep(cols['G'],currow['G']),
                 rep(cols['B'],currow['B']))
        r2s <- sample(r2s)
        bc1 <- substr(bc,1,8)
        bc2 <- substr(bc,9,14)
        umis <- replicate(nUMI*2,{ paste0(sample(c('A','C','G','T'),9,replace=TRUE),collapse='')})
        umis <- unique(umis)[1:nUMI]
        r1s <- paste0(bc1,linker,bc2,umis,'T',sep='')
        
        dels <- replicate(nUMI,{ sample(c(TRUE,FALSE),42,prob=c((1-dr),dr),replace=TRUE) })
        spl <- strsplit(r1s,'')
        r1s_wdels <- sapply(1:nUMI,\(i){
          x <- c(spl[[i]][dels[,i]],ts)[1:42]
          return(paste0(x,collapse=''))
        })
        
        id <- paste0(baseid,3000+(100*k)+1:nUMI,'/')
        if (i==1){
          write(paste0(id,'1\n',r1s,'\n+\n',r1qual,collapse='\n'),paste0('./data/',imgName,'_',k,'_r1_perfect.fq'),append=FALSE)
          write(paste0(id,'1\n',r1s_wdels,'\n+\n',r1qual,collapse='\n'),paste0('./data/',imgName,'_',k,'_r1.fq'),append=FALSE)
          write(paste0(id,'2\n',r2s,'\n+\n',r2qual,collapse='\n'),paste0('./data/',imgName,'_',k,'_r2.fq'),append=FALSE)
        } else {
          write(paste0(id,'1\n',r1s,'\n+\n',r1qual,collapse='\n'),paste0('./data/',imgName,'_',k,'_r1_perfect.fq'),append=TRUE)
          write(paste0(id,'1\n',r1s_wdels,'\n+\n',r1qual,collapse='\n'),paste0('./data/',imgName,'_',k,'_r1.fq'),append=TRUE)
          write(paste0(id,'2\n',r2s,'\n+\n',r2qual,collapse='\n'),paste0('./data/',imgName,'_',k,'_r2.fq'),append=TRUE)
        }
      }
    }
    system(paste0('cat ./data/',imgName,'_*_r1_perfect.fq > ./data/',imgName,'_r1_perfect.fq'))
    system(paste0('cat ./data/',imgName,'_*_r1.fq > ./data/',imgName,'_r1.fq'))
    system(paste0('cat ./data/',imgName,'_*_r2.fq > ./data/',imgName,'_r2.fq'))
    system(paste0('rm ./data/',imgName,'_*_r*.fq'))
  }
}

# Now the Read 1 FASTQs can be processed with the
# Syrah package and both versions run through the
# WARP Slide-seq pipeline using the default 
# (workspace) mouse reference. Use mmus_seqs.csv
# to translate gene IDs to R/G/B and the 
# *_coords.tsv file for bead barcodes/coordinates.