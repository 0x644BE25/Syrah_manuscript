######################################################
# GENERATE PLOTS AND XLSX
#
# GOAL: Summarize data and generate all plots.
######################################################

setwd('/n/projects/cb2350/syrah/v2/code/manuscript/')

# ================= IMPORTS ==========================

pkgs_needed <- setdiff(c('ggplot2','ggrepel','patchwork','openxlsx'),installed.packages())
install.packages(pkgs_needed)

library(ggplot2)
library(ggrepel)
library(patchwork)
library(openxlsx)

# ================= PARAMS ===========================

datasets <- c('chick','planarian 1 run','planarian 2 runs','Curio test data')
versions <- c('standard','Syrah')
measures <- c('genes','beads','UMIs')

dataset_colors <- setNames(c('#77cc33','#8aacdb','#638bc7','#ee8866'),datasets)

planarian_colors <- c('base'='#f1c267','additional sequencing'='#e3a13e',
                      'Syrah'='#8aacdb','both'='#638bc7')
read_cols <- c('FALSE'='#bbbbbb','TRUE'='#888888',
               '1'='white','2'='burlywood','3'='darkorange','4'='yellow',
               '5'='chartreuse','6'='forestgreen','7'='cyan2','8'='dodgerblue',
               '9'='darkorchid','10'='hotpink1','11'='brown3','12'='black')

eb <- element_blank()

ptsize <- 0.1353801

# ================= COLLATE DATA =====================

files <- dir('./data/')

df <- do.call(rbind,lapply(df,function(f){
  counts <- read.delim(paste0('./data/',f),row.names=1,header=TRUE)
  
  beads <- ncol(counts)
  genes <- sum(rowSums(counts)>0)
  UMIs <- sum(colSums(counts))
  
  version <- tail(head(strsplit(f,'\\_')[[1]],-2),1)
  
  return(data.frame(dataset,version,beads,genes,UMIs))
}))
df$dataset <- factor(df$dataset,levels=datasets)
df$version <- factor(df$version,levels=versions)
write.csv(df,'./data/temp.csv',row.names=1)

# ================= FIG3A: COMPARE DATASETS ==========

df3a <- do.call(rbind,lapply(datasets,function(ds){
  syrah <- df[df$version=='Syrah' & df$dataset==ds,]
  standard <- df[df$version=='standard' & df$dataset==ds,]
  
  # fold change
  syrah[,measures] <- syrah[,measures]/standard[,measures]
  
  # pct change
  syrah[,measures] <- 100*(syrah[,measures]-1)
  
  return(syrah[,c('dataset',measures)])
}))

fig3a <- ggplot(reshape2::melt(df3a),aes(x=dataset,fill=dataset,y=value)) + 
  geom_col() + 
  scale_fill_manual(values=dataset_colors,labels=\(x){sub(' ','\n',x)},name='') +
  ylab('% change vs. standard') + 
  theme_bw() + 
  theme(axis.title.x=eb,axis.text.x=eb,axis.ticks.x=eb,panel.grid.major.x=eb,
        panel.grid.major.y=element_line(color='#CCCCCC',linetype='dashed',linewidth=1),
        panel.grid.minor=eb,
        legend.position='right',legend.key.spacing.y=unit(.1,'in'),legend.title=eb) +
  guides(fill=guide_legend(byrow=TRUE)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,30))  +
  facet_grid(.~factor(variable,levels=measures))
ggsave(paste0('./fig3a_compare_datasets.pdf'),fig3a,width=(8*.75),height=(2.5*.75),units='in',dpi=100)

# ================= FIG1B: SYRAH VS STANDARD =========

std1run <- df[df$dataset=='planarian 1 run' & df$version=='standard',measures]
addlseq <- df[df$dataset=='planarian 2 runs' & df$version=='standard',measures]
syrah <- df[df$dataset=='planarian 1 run' & df$version=='Syrah',measures]
both <- df[df$dataset=='planarian 2 runs' & df$version=='Syrah',measures]

# compare to standard 1 run
df3b <- rbind(
  data.frame(type='additional sequencing',(100*((addlseq/std1run)-1))),
  data.frame(type='Syrah',(100*((syrah/std1run)-1))),
  data.frame(type='both',(100*((both/std1run)-1)))
)
df3b$type <- factor(df3b$type,levels=df3b$type)

fig3b <- ggplot(reshape2::melt(df3b),aes(x=type,fill=type,y=value)) +
  geom_col() +
  scale_fill_manual(values=planarian_colors,name='',labels=\(x){ gsub(' ','\n',x) }) + 
  facet_grid(.~factor(variable,levels=measures)) +
  theme_bw() + 
  theme(axis.title.x=eb,axis.text.x=eb,axis.ticks.x=eb,panel.grid.major.x=eb,
        panel.grid.major.y=element_line(color='#CCCCCC',linetype='dashed',linewidth=1),
        panel.grid.minor=eb,
        legend.position='right',strip.text=element_text(size=rel(1)),
        legend.key.spacing.y=unit(.1,'in'),legend.title=eb) +
  guides(fill=guide_legend(byrow=TRUE)) +
  scale_y_continuous(expand=expansion(mult=c(0,.05)),breaks=c(10,30,50)) +
  labs(x='',y='% increase')

ggsave('./fig3b_Syrah_vs_sequncing.pdf',fig3b,width=(8*.75),height=(2.5*.75),units='in',dpi=100)

# ================= FIGS1A: EXAMPLE READS ============

beads_standard <- as.character(read.delim('./data/planarian_2_runs_standard_counts_min10umi.tsv.gz',nrows=1,header=FALSE))
beads_syrah <- as.character(read.delim('./data/planarian_2_runs_Syrah_counts_min10umi.tsv.gz',nrows=1,header=FALSE))

puck <- read.delim('./data/planarian_bead_coordinates.txt',sep='\t',header=TRUE,row.names=1); colnames(puck) <- c('x','y')
puck <- puck[unique(c(beads_standard,beads_syrah)),]

# under tissue beads manually selected with custom ShinyApp
puck$under_tissue <- read.csv('./data/planarian_bead_under_tissue_designations.csv',row.names=1)[rownames(puck),'under_tissue']
puck$standard <- rownames(puck) %in% beads_standard
puck$syrah <- rownames(puck) %in% beads_syrah

# example reads and their beads
reads <- read.csv('./data/planarian_example_reads.csv')
reads$n <- factor(reads$n,levels=1:12)

# standard version
pstd <- ggplot(data=NULL,aes()) + geom_point(data=puck[puck$standard,],aes(x=x,y=y,color=under_tissue),size=ptsize) +
  scale_color_manual(values=read_cols) +
  geom_label_repel(data=reads,aes(x=orig.x,y=orig.y,fill=n,label=n),color='black',min.segment.length=0) +
  geom_point(data=reads,aes(x=orig.x,y=orig.y,fill=n),size=1,pch=21) + 
  scale_fill_manual(values=read_cols) +
  theme_void() + coord_fixed() + 
  theme(legend.position='none')

# Syrah version
psyr <- ggplot(data=NULL,aes()) + geom_point(data=puck[puck$syrah,],aes(x=x,y=y,color=under_tissue),size=ptsize) +
  scale_color_manual(values=read_cols) +
  geom_label_repel(data=reads,aes(x=syrah.x,y=syrah.y,fill=n,label=n),color='black',min.segment.length=0) +
  geom_point(data=reads,aes(x=syrah.x,y=syrah.y,fill=n),size=1,pch=21) + 
  scale_fill_manual(values=read_cols) +
  theme_void() + coord_fixed() + 
  theme(legend.position='none')

pS1a <- patchwork::wrap_plots(list(pstd,psyr),ncol=2)
ggsave('./figS1a_example_reads.png',pS1a,width=7,height=3.5,units='in',dpi=300)

# ================= DATA XLSX ========================

# more informative column names
colnames(df3a) <- c('dataset',paste0('pct_change_',measures))
colnames(df3b) <- c('type',paste0('pct_change_',measures))

xlsx <- list('dataset_stats'=df,'fig3a'=df3a,'fig3b'=df3b,'figS1a'=reads)
openxlsx::write.xlsx(xlsx,'./figure_data.xlsx')
