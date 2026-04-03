# THIS IS TO RUN SYRAHS COMPONENTS MULTIPLE TIMES FOR BENCHMARKING PURPOSES
if (!require('Syrah')) {
  if (!require('devtools')) { install.packages('devtools') }
  devtools::install_github('0x644BE25/Syrah',force=TRUE)
}
if (!require('peakRAM')) {
  if (!require('devtools')) { install.packages('devtools') }
  devtools::install_github('tpq/peakRAM')
}
library('Syrah')
library('peakRAM')

# PARAMS
n <- 5
outFile <- 'benchmarking_results.csv'
datasets <- list('Curio'=c(url='https://github.com/0x644BE25/Syrah/raw/refs/heads/main/test_data/example_input_mouse_spleen_1M.tar.gz',
                           bcs='A0010_039_BeadBarcodes.txt',
                           r1fq='Mouse_spleen_1m_R1.fastq.gz'),
                 'SlideSeq'=c(url='https://github.com/0x644BE25/Syrah/raw/refs/heads/main/test_data/planarian_test_data.tar.gz',
                                 bcs='planarian_L43430_coordinates.txt',
                                 r1fq='planarian_L43430_H5F73AFX2_1M_R1.fastq.gz'))

# GET TEST DATA
for (ds in names(datasets)) {
  download.file(datasets[[ds]]['url'],paste0(ds,'.tar.gz'))
  untar(paste0(ds,'.tar.gz'))
}

# RUN SYRAH
res <- data.frame()
cat('\nSyrah benchmarking begins!')
for (ds in names(datasets)) {
  bcs <- datasets[[ds]]['bcs']
  r1fq <- datasets[[ds]]['r1fq']
  dd <- paste0(bcs,'_dedup_map.txt')
  wl <- paste0(bcs,'_whitelist.txt')
  for (i in 1:n) {
    cat('\n ===================',ds,' rep',i,'/',n,'============================\n')
    a <- peakRAM(make_bead_dedup_map(bcs))
    b <- peakRAM(make_barcode_whitelist(dd,bcs))
    c <- peakRAM(correct_barcodes(wl,r1fq))
    curr <- rbind(a,b,c)
    curr$rep <- i
    curr$datset <- ds
    res <- rbind(res,curr)
    write.csv(res,outFile,row.names=FALSE)
  }
}
cat('\n\nALL DONE! Thanks for your help :)\n')
cat('I just need the',outFile,'file, please. Everything else can be deleted.\n\n')
