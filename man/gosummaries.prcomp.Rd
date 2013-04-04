\name{gosummaries.prcomp}
\alias{gosummaries.prcomp}
\title{Prepare gosummaries object based on PCA results}
\usage{
  \method{gosummaries}{prcomp} (x, annotation = NULL,
    components = 1:6, n_genes = 500, ...)
}
\arguments{
  \item{x}{an object of class \code{prcomp}}

  \item{annotation}{a \code{data.frame} describing the
  samples, its row names should match with column names of
  the projection matrix in x}

  \item{components}{numeric vector of components to
  include}

  \item{n_genes}{how many top genes to use for annotation}

  \item{\dots}{GO annotation filtering parameters as
  defined in \code{\link{gosummaries.default}}}
}
\value{
  A gosummaries object.
}
\description{
  The PCA results are converted into a gosummaries object,
  by extracting genes with the largest positive and
  negative weights from each component.
}
\details{
  The usual visualisation of PCA results displays the
  projections of sample expression on the principal axes.
  It shows if and how the samples cluster, but not why do
  they behave like that. Actually, it is possible to go
  further and annotate the axes by studying genes that have
  the largest influence in the linear combinations that
  define the principal components. For example, high
  expression of genes with large negative weights pushes
  the samples projection to the negative side of the
  principal axis and large positive weigths to the positive
  side. If a sample has highly expressed genes in both
  groups it stays most probably in the middle. If we
  annotate functionally the genes with highest positive and
  negative weights for each of the principal axes, then it
  is possible to say which biological processes drive the
  separation of samples on them.

  This function creates a gosummaries object for such
  analysis. It expects the results of \code{\link{prcomp}}
  function. It assumes that the PCA was done on samples
  and, thus, the row names of the rotation matrix can be
  interpreted as gene names. For each component it
  annotates \code{n_genes} elements with highest positive
  and negative weights.
}
\examples{
\dontrun{
data(tissue_example)

pcr = prcomp(t(tissue_example$exp))
gs_pca = gosummaries(pcr, annotation = tissue_example$annot)

plot(gs_pca, classes = "Tissue")
}
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

