These are the (numpy) sparse matrix components from the WARP pipeline.
Even with sparse matrix representation, these are much smaller than the matrix itself.
You can easily re-constitute them to sparse matrices in R with the getMatrixFromNpz 
function from the helper_functions.R script in the main repo directory.

Human is split into pieces due to size constraints, but all other files
include both original and Syrah samples as well as the coordinates file.

The planarian data also includes pre-computed cell type prediction (LabelTransfer) data,
as it requires extensive preparation of data from Benham-Pyle et al. 2021.
See additional_sequencing_DE_analysis_and_cell_type_predictions.R for more information 
on performing this analysis yourself.
