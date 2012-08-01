\name{panel_boxplot}
\alias{panel_boxplot}
\alias{panel_violin}
\alias{panel_violin_box}
\title{Panel drawing functions}
\usage{
  panel_boxplot(data, fontsize = 10, par)

  panel_violin(data, fontsize = 10, par)

  panel_violin_box(data, fontsize = 10, par)
}
\arguments{
  \item{data}{the data from Data slot of the go sumaries
  object}

  \item{fontsize}{fontsize in points, mainly used to ensure
  that the legend fontsizes match}

  \item{par}{additional parameters for drawing the plot,
  given as list. These functions use only \code{classes}
  slot for figuring out which parameter to use for coloring
  the "geom"-s. However, when building a custom function it
  provides a way to give extra parameters to the plotting
  function.}
}
\value{
  It returns a function that can draw a ggplot2 plot of the
  data in Data slot of a gosummaries object. The legend and
  the actual plots for the panels are extracted later from
  the figure produced by this function.
}
\description{
  These functions are used to draw the panel portion of
  every component based on the Data slots in gosummaries
  object. These concrete functions assume the data is
  presented as is done by
  \code{\link{add_expression.gosummaries}}. They provide
  three options: boxplot, violin plot (which shows the
  distrubution more precisely) and both combined.
}
\details{
  These functions specify in principle the general setting
  for the figures, like which "geom"-s, how the data is
  transformed and summarized, etc. To make small
  adjustments to the figure such as changeing color scheme,
  write your own customization function (See
  \code{\link{customize}} as example).

  It is possible to write your own panel plotting function,
  as long as the parameters used and the return value are
  similar to what is specified here. When writing a new
  panel function one only has to make sure that it matches
  the data given in the Data slot of the gosummaries
  object.
}
\examples{
# Example1
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

