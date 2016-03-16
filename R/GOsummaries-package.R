#' Word cloud summaries of GO enrichment analysis
#'
#' A package to visualize Gene Ontology (GO) enrichment analysis results on gene lists 
#' arising from different analyses such clustering or PCA. The significant GO categories 
#' are visualised as word clouds that can be combined with different plots summarizing 
#' the underlying data.  
#' 
#' The goal of GOsummaries package is to draw figures that can be used in 
#' presentations and articles. To draw them, the user should first construct a 
#' \code{\link{gosummaries}} object and then use its plot function on it. One can start 
#' constructing the \code{\link{gosummaries}} object from gene lists, with filling in all 
#' the necessary information step by step. However, there are some convenience functions 
#' for different classes of common analysis results. See \code{\link{gosummaries.kmeans}}, 
#' \code{\link{gosummaries.MArrayLM}} and \code{\link{gosummaries.prcomp}} corresponding 
#' to k-means, limma and PCA results.
#' 
#' The \code{\link{plot.gosummaries}} describes how to customize the plots. 
#' 
#' The word cloud drawing function \code{\link{plotWordcloud}} in this package is 
#' implemented largely based on the code from package \code{wordcloud}, but with slight 
#' tweaks: it uses \code{grid} graphics, has some additional layout options, has more 
#' intelligent options to scale the text sizes to fit the picture and, finally, should be 
#' a bit faster since larger part of the algorithm was implemeted in C++. 
#' 
#' 
#' @name GOsummaries-package
#' @docType package
#' 
#' @import grid
#' @import plyr
#' @import gProfileR
#' @import reshape2
#' @import limma
#' @import gtable
#' @import Rcpp
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggplot2 ggplot_gtable
#' @importFrom ggplot2 layer
#' @importFrom ggplot2 qplot
#' @importFrom ggplot2 scale_fill_discrete
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 position_identity
#' @importFrom ggplot2 scale_fill_brewer
#' 
#' @useDynLib GOsummaries 
#' 
#' @aliases GOsummaries
NULL