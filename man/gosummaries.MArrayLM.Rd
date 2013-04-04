\name{gosummaries.MArrayLM}
\alias{gosummaries.MArrayLM}
\title{Prepare gosummaries object based on limma results}
\usage{
  \method{gosummaries}{MArrayLM} (x, p.value = 0.05,
    lfc = 1, adjust.method = "fdr", exp = NULL,
    annotation = NULL, components = 1:ncol(x), ...)
}
\arguments{
  \item{x}{an object of class \code{MArrayLM}}

  \item{p.value}{p-value threshold as defined in topTable}

  \item{lfc}{log fold change threshold as defined in
  topTable}

  \item{adjust.method}{multiple testing adjustment method
  as defined in topTable}

  \item{exp}{an expression matrix, with row names
  corresponding to the names of the genes in clusters
  (Optional)}

  \item{annotation}{a \code{data.frame} describing the
  samples, its row names should match with column names of
  \code{exp} (Optional)}

  \item{components}{numeric vector of comparisons to
  annotate}

  \item{\dots}{GO annotation filtering parameters as
  defined in \code{\link{gosummaries.default}}}
}
\value{
  A gosummaries object.
}
\description{
  The gosummaries object is created based on the
  differentially expresed genes, each contrast defines one
  component.
}
\details{
  The usual differential expression analysis involves
  making several comparisons between treatments ehere each
  one yields an up and down regulated gene list. In a
  GOsummaries figure each comparison is displayed as one
  component with two wordclouds. If expression matrix is
  attached then the panel shows the expression values for
  each gene as boxplots, if not then number of genes is
  displayed
}
\examples{
\dontrun{
data(tissue_example)

# Do the t-test comparisons
mm = model.matrix(~ factor(tissue_example$annot$Tissue) - 1)
colnames(mm) = make.names(levels(factor(tissue_example$annot$Tissue)))

contrast = limma::makeContrasts(brain - cell.line, hematopoietic.system - muscle, cell.line - hematopoietic.system, levels = colnames(mm))

fit = limma::lmFit(tissue_example$exp, mm)
fit = limma::contrasts.fit(fit, contrast)
fit = limma::eBayes(fit)

gs_limma = gosummaries(fit)
gs_limma_exp = gosummaries(fit, exp = tissue_example$exp, annotation = tissue_example$annot)

plot(gs_limma, fontsize = 8)
plot(gs_limma, panel_height = 0, fontsize = 8)
plot(gs_limma_exp, classes = "Tissue", fontsize = 8)
}
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

