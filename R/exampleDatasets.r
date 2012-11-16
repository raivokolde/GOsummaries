#' \code{tissue_example} is a dataset based on Lukk \emph{et al}, it contains a subset of 24 samples 
#' from more than 5000 in the original article. The gene expression in all the samples is measured using 
#' Affymetrix U133A array. 
#' The dataset is a list with 2 slots
#' \enumerate{ 
#' 	\item \code{exp} - expression matrix;
#' 	\item \code{annot} - annotation data frame.
#' }
#' 
#' The GOsummaries objects based on this data created in different examples are alaso attached to the 
#' package. These can be used to test the package if no internet connection is available. Names of these are \code{gs_kmeans}, \code{gs_limma}, \code{gs_limma_exp} and \code{gs_pca} 
#' 
#' @name tissue_example 
#' 
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' 
#' @references Lukk M, Kapushesky M, Nikkila J, Parkinson H, Goncalves A, Huber W, 
#' Ukkonen E, Brazma A. "A global map of human gene expression." Nat Biotechnology. 2010 
#' Apr;28(4):322-4.
#' 
#' @aliases gs_kmeans gs_limma gs_limma_exp gs_pca
#' 
#' @keywords data
NULL
