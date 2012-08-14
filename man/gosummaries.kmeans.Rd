\name{gosummaries.kmeans}
\alias{gosummaries.kmeans}
\title{Prepare gosummaries object based on k-means results}
\usage{
  \method{gosummaries}{kmeans} (x, exp = NULL,
    annotation = NULL, components = 1:length(x$size), ...)
}
\arguments{
  \item{x}{an object of class \code{kmeans}}

  \item{exp}{an expression matrix, with row names
  corresponding to the names of the genes in clusters
  (Optional)}

  \item{annotation}{a \code{data.frame} describing the
  samples, its row names should match with column names of
  \code{exp} (Optional)}

  \item{components}{numeric vector of clusters to annotate}

  \item{\dots}{GO annotation filtering parameters as
  defined in \code{\link{gosummaries.default}}}
}
\value{
  A gosummaries object.
}
\description{
  The gosummaries object is created based on the genes in
  the clusters, it is possible to add corresponding gene
  expression data as well.
}
\details{
  The k-means clustering of expression matrix naturally
  defines a set of gene lists that can be annotated
  functionally and displayed as a GOsummaries figure. This
  functon takes in a \code{kmeans} object and and converts
  it to a \code{gosummaries} object that can be plotted. If
  expression matrix is attached then the panel shows the
  expression values for each gene as boxplots, if not then
  number of genes is displayed

  It is advisable to filter some genes out before doing the
  clustering since the very large gene lists (more than
  2000 genes) might fail the annotation step and are
  usually not too specific either.
}
\examples{
# Example1
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

