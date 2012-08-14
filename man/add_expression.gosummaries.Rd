\name{add_expression.gosummaries}
\alias{add_expression.gosummaries}
\title{Add expression data to gosummaries object}
\usage{
  add_expression.gosummaries(gosummaries, exp, annotation)
}
\arguments{
  \item{gosummaries}{a gosummaries object}

  \item{exp}{an expression matrix, with row names
  corresponding to the names in the Gene_lists slot}

  \item{annotation}{a \code{data.frame} describing the
  samples, its row names should match with column names of
  \code{exp}}
}
\description{
  Function to add expression data and its annotations to a
  gosummaries object.
}
\details{
  The data is added to the object in a "long" format so it
  would be directly usable by the ggplot2 based panel
  drawing functions \code{\link{panel_boxplot}} etc. For
  each component it produces a data frame with columns:
  \itemize{ \item \bold{x} : sample IDs for the x axis,
  their factor order is the same as on the columns of input
  matrix \item \bold{y} : expression values from the matrix
  \item . . . : sample annotation columns from the
  annotation table that can be displayed on figure as
  colors. }
}
\examples{
# Example1
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

