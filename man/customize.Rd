\name{customize}
\alias{customize}
\title{Customization function for panel}
\usage{
  customize(p, par)
}
\arguments{
  \item{p}{a ggplot2 plot object}

  \item{par}{parameters object like in
  \code{\link{panel_boxplot}}}
}
\value{
  a ggplot2 plot object with added customizations
}
\description{
  This function is supposed to make small changes in the
  panel function appearance like changing color scheme for
  example. It has to match with the output of the
  corresponding panel function. Check examples in
  plot.gosummaries to see how to write one yourself.
}
\examples{
# Example1
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

