######################################################
# HELPER FUNCTIONS
######################################################


if (!require('callr')) install.packages('callr')
if (!require('reticulate')) install.packages('reticulate')

getCountsFromNpz <- function(f) {
  library('callr')
  counts <- callr::r(function(x) {
    library('reticulate')
    library('Matrix')
    reticulate::use_condaenv('/home/cb2350/miniforge3')
    np <- reticulate::import('numpy')
    sp <- reticulate::import('scipy.sparse')
    counts <- as(sp$load_npz(x),"CsparseMatrix")
    colnames(counts) <- as.character(np$load(gsub('.npz','_col_index.npy',x)))
    rownames(counts) <- as.character(np$load(gsub('.npz','_row_index.npy',x)))
    return(Matrix::t(counts))
  },args=list(f))
  return(counts)
}
