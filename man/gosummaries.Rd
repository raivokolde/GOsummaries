\name{gosummaries}
\alias{gosummaries}
\alias{gosummaries.default}
\title{Constructor for gosummaries object}
\usage{
  gosummaries(x, ...)

  \method{gosummaries}{default} (x, organism = "hsapiens",
    go_branches = c("BP", "KEGG", "REAC"),
    p_value_threshold = 0.01, min_size = 50,
    max_size = 1000, max_signif = 40, ...)
}
\arguments{
  \item{x}{list of arrays of gene names (or list of lists
  of arrays of gene names)}

  \item{organism}{the organism that the gene lists
  correspond to. The format should be as follows:
  "hsapiens", "mmusculus", "scerevisiae", etc.}

  \item{\dots}{Description of parameter}

  \item{go_branches}{GO tree branches and pathway databases
  as denoted in g:Profiler (Possible values: BP, CC, MF,
  KEGG, REAC)}

  \item{p_value_threshold}{threshold for p-values that are
  corrected for multiple testing}

  \item{min_size}{minimal size of functional category to be
  considered}

  \item{max_size}{maximal size of functional category to be
  considered}

  \item{max_signif}{maximal number of categories returned
  per query}
}
\value{
  A gosummaries type of object
}
\description{
  Constructor for gosummaries object that contains all the
  necessary information to d raw the figure, like gene
  lists and their annotations, expression data and all the
  relevant texts.
}
\details{
  The object is a list of "components", with each component
  defined by a gene list or a pair of gene lists. Each
  "component" has the slots as follows: \itemize{ \item
  \bold{Title}: title string of the component. (Default:
  the names of the gene lists) \item \bold{Gene_lists}:
  list of one or two gene lists with names gl1 (and gl2 if
  present).  \item \bold{GPR}: g:Profiler results based on
  the Gene_lists slot.  \item \bold{Data}: the related data
  (expression values, PCA rotation, ...) that is used to
  draw the "panel" i.e. theplot above the wordclouds. In
  principle there is no limitation what kind of data is
  there as far as the function that is provided to draw
  that in \code{\link{plot.gosummaries}} can use it.  \item
  \bold{Percentage}: a text that is drawn on the right top
  corner of every component. In case of PCA this is the
  percentage of variation the component explains, by
  default it just depicts the number of genes in the
  Gene_lists slot. }

  The GO enrichment analysis is performed using g:Profiler
  web toolkit and its associated R package
  \code{gProfileR}. This means the computer has to have
  internet access to annotate the gene lists. Since
  g:Profiler can accept a wide range of gene IDs then user
  usually does not have to worry about converitng the gene
  IDs into right format. To be absolutely sure the tool
  recognizes the gene IDs one can check if they will give
  any results in
  \url{http://biit.cs.ut.ee/gprofiler/gconvert.cgi}.

  There can be a lot of results for a typical GO enrichment
  analysis but usually these tend to be pretty redundant.
  Since one can fit only a small number of categories into
  a word cloud we have to bring down the number of
  categories to show an reduce the redundancy. For this we
  apply an algortihm that selects from groups of related
  categories only the ones with the best p-values. The
  algorithm works as follows. \itemize{ \item Delete all
  the results where category size is smaller than
  \code{min_size}, larger than \code{max_size}, enrichment
  p-value is larger than \code{p_value_threshold} or the
  category does not belong to one of the
  \code{go_branches}.  \item Map all the remaining results
  back to the GO graph and find connected components.
  \item From each connected component retain the category
  with smallest p-value.  \item If more than
  \code{max_signif} categories are still present return
  \code{max_signif} ones with the best p-values }

  The slots of the object can be filled with custom
  information using a function
  \code{\link{add_to_slot.gosummaries}}.

  By default the Data slot is filled with a dataset that
  contains the number of genes in the Gene_lists slot.
  Expression data can be added to the object for example by
  using function \code{\link{add_expression.gosummaries}}.
  It is possible to derive your own format for the Data
  slot as well, as long as a panel plotting function for
  this data is alaso provided (See
  \code{\link{panel_boxplot}} for further information).

  There are several constructors of gosummaries object that
  work on common analysis result objects, such as
  \code{\link{gosummaries.kmeans}},
  \code{\link{gosummaries.MArrayLM}} and
  \code{\link{gosummaries.prcomp}} corresponding to
  k-means, limma and PCA results.
}
\examples{
# Example1
# Example1
}
\author{
  Raivo Kolde <rkolde@gmail.com>

  Raivo Kolde <rkolde@gmail.com>
}
\seealso{
  \code{\link{gosummaries.kmeans}},
  \code{\link{gosummaries.MArrayLM}},
  \code{\link{gosummaries.prcomp}}
}

