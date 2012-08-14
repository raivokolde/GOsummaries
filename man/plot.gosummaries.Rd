\name{plot.gosummaries}
\alias{plot.gosummaries}
\title{Plot the GOsummaries figure}
\usage{
  \method{plot}{gosummaries} (x,
    components = names(x)[1:min(10, length(x))],
    panel_plot = NULL, panel_customize = NULL,
    panel_par = list(), panel_height = 5, fontsize = 10,
    term_length = 35,
    wordcloud_colors = c("grey70", "grey10"),
    filename = NA, ...)
}
\arguments{
  \item{x}{a gosummaries object}

  \item{components}{names of the gosummaries components to
  draw. The names are usually numbers 1:n in character
  format.}

  \item{panel_plot}{plotting function for panel}

  \item{panel_customize}{customization function for the
  panel plot, menat for making small changes like changing
  colour scheme}

  \item{panel_par}{list of arguments passed on to
  \code{panel_plot} function}

  \item{panel_height}{panel height as number of lines, with
  given \code{fontsize}. If set to 0 no panel is drawn.}

  \item{fontsize}{font size used throughout the figure in
  points}

  \item{term_length}{maximum length of the dispalyed GO
  categories in characters, longer names are cropped to
  this size}

  \item{wordcloud_colors}{two element vector of colors to
  define color scheme for displaying the enrichment
  p-values across the wordclouds. First element defines the
  color for category with worst p-value and the second for
  the word with the best.}

  \item{filename}{file path where to save the picture.
  Filetype is decided by the extension in the path.
  Currently following formats are supported: png, pdf,
  tiff, bmp, jpeg. Even if the plot does not fit into the
  plotting window, the file size is calculated so that the
  plot would fit there.}

  \item{\dots}{not used}
}
\description{
  The function to draw a GOsummaries figure based on a
  \code{gosummaries} object.  The GOsummaries figure
  consists of several components each defined by a gene
  list ora a pair of them. The GO annotations of them are
  shown as wordclouds. Optionally one can draw related
  (expression) data on panels atop of the wordclouds.
}
\details{
  In most cases the function can decide which type of plot
  to draw into the panel part. If there is no data
  explicitly put into the Data slots of the gosummaries
  object, it just draws a horizontal barplot with the
  numbers of genes. On visualizing the PCA data it draws
  histogram of the samples on the principal axes. For
  clustering and differential expression it draws the
  boxplot of expression values.
}
\examples{
# Example1
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

