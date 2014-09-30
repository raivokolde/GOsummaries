#' Example gene expression dataset 
#' 
#' \code{tissue_example} is a dataset extracted from Lukk \emph{et al}, it contains a subset of 24 samples 
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
#' @docType data
#' 
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' 
#' @references Lukk M, Kapushesky M, Nikkila J, Parkinson H, Goncalves A, Huber W, 
#' Ukkonen E, Brazma A. "A global map of human gene expression." Nat Biotechnology. 2010 
#' Apr;28(4):322-4.
#' 
#' @aliases gs_kmeans gs_limma gs_limma_exp gs_pca
#' 
#' @keywords data
NULL


#' Example metabolomic dataset
#' 
#' \code{metabolomic_example} is a dataset extracted from York \emph{et al} (Metabolights ID:
#'  MTBLS30), it contains a subset of 120 wild-type samples from 4 tissues: heart, skeletal
#'  muscle, liver and brain.
#' 
#' The dataset is a list with 2 slots
#' \enumerate{ 
#' 	\item \code{data} - metabolite concentration matrix;
#' 	\item \code{annot} - annotation data frame.
#' }
#' 
#' @name metabolomic_example 
#' @docType data
#' 
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' 
#' @references York, B., Sagen, J. V., Tsimelzon, A., Louet, J. F., Chopra, A. R., Reineke, E. L.,
#'  et al. (2013). Research resource: tissue- and pathway-specific metabolomic profiles of the
#'  steroid receptor coactivator (SRC) family. Mol Endocrinol, 27(2), 366-380.
#' 
#' @keywords data
NULL

#' Example metagenomic dataset
#' 
#' \code{metagenomic_example} is a small sample of Human Microbiome Project 16S dataset for
#'  finding biomarkers characterizing different level of oxygen availability in different
#'  bodysites. This was downloaded from http://huttenhower.sph.harvard.edu/webfm_send/129. 
#' 
#' The dataset is a list with 2 slots
#' \enumerate{ 
#' 	\item \code{otu} - otu table for this experiment;
#' 	\item \code{annot} - annotation data frame.
#' }
#' 
#' @name metagenomic_example 
#' @docType data
#' 
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' 
#' @keywords data
NULL