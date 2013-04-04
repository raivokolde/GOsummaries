\name{gosummaries}
\alias{gosummaries}
\alias{gosummaries.default}
\title{Constructor for gosummaries object}
\usage{
  gosummaries(x, ...)

  \method{gosummaries}{default} (x, organism = "hsapiens",
    go_branches = c("BP", "ke", "re"), max_p_value = 0.01,
    min_set_size = 50, max_set_size = 1000,
    max_signif = 40, ordered_query = T,
    hier_filtering = "moderate", ...)
}
\arguments{
  \item{x}{list of arrays of gene names (or list of lists
  of arrays of gene names)}

  \item{organism}{the organism that the gene lists
  correspond to. The format should be as follows:
  "hsapiens", "mmusculus", "scerevisiae", etc.}

  \item{\dots}{additional parameters for gprofiler
  function}

  \item{go_branches}{GO tree branches and pathway databases
  as denoted in g:Profiler (Possible values: BP, CC, MF,
  ke, re)}

  \item{max_p_value}{threshold for p-values that are
  corrected for multiple testing}

  \item{min_set_size}{minimal size of functional category
  to be considered}

  \item{max_set_size}{maximal size of functional category
  to be considered}

  \item{max_signif}{maximal number of categories returned
  per query}

  \item{ordered_query}{logical showing if the lists are
  ordered or not (it determines if the ordered query
  algorithm is used in g:Profiler)}

  \item{hier_filtering}{a type of hierarchical filtering
  used when reducing the number of g:Profiler results (see
  \code{\link{gprofiler}} for further information)}
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
  use hierarchical filtering option \"moderate\" in
  g:Profiler. In g:Profiler the categories are grouped
  together when they share one or more enriched parents.
  The \"moderate\" option selects the most significant
  category from each of such groups. (See more at
  http://biit.cs.ut.ee/gprofiler/)

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
\dontrun{
# Define gene lists
genes1 = c("203485_at", "209469_at", "209470_s_at", "203999_at", "205358_at", "203130_s_at",
"210222_s_at", "202508_s_at", "203001_s_at", "207957_s_at", "203540_at", "203000_at", "219619_at",
"221805_at", "214046_at", "213135_at", "203889_at", "209990_s_at", "210016_at", "202507_s_at",
"209839_at", "204953_at", "209167_at", "209685_s_at",  "211276_at", "202391_at", "205591_at",
"201313_at")
genes2 = c("201890_at", "202503_s_at", "204170_s_at", "201291_s_at", "202589_at", "218499_at",
"209773_s_at", "204026_s_at", "216237_s_at", "202546_at", "218883_s_at", "204285_s_at",
"208659_at", "201292_at", "201664_at")


gl1 = list(List1 = genes1,  List2 = genes2) # One list per component
gl2 = list(List = list(genes1, genes2)) # Two lists per component

# Construct gosummaries objects
gs1 = gosummaries(gl1)
gs2 = gosummaries(gl2)

plot(gs1, fontsize = 8)
plot(gs2, fontsize = 8)

# Changing slot contents using using addToSlot.gosummaries
gs1 = add_to_slot.gosummaries(gs1, "Title", list("Neurons", "Cell lines"))

# Adding expression data
data(tissue_example)

gs1 = add_expression.gosummaries(gs1, exp = tissue_example$exp, annotation = tissue_example$annot)
gs2 = add_expression.gosummaries(gs2, exp = tissue_example$exp, annotation = tissue_example$annot)

plot(gs1, panel_par = list(classes = "Tissue"), fontsize = 8)
plot(gs2, panel_par = list(classes = "Tissue"), fontsize = 8)
}
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

