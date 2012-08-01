\name{plotWordcloud}
\alias{plotWordcloud}
\title{Plot a wordcloud given words and frequencies}
\usage{
  plotWordcloud(words, freq, rot.per = 0.3,
    max_min = c(1, 0.1), scale = 0.4, min.freq = 3,
    max.words = Inf, random.order = FALSE,
    colors = "black", random.colors = FALSE,
    algorithm = "circle", tryfit = T, add = F)
}
\arguments{
  \item{words}{vector of words to draw}

  \item{freq}{frequencies for words, has to be the same
  length as words vector}

  \item{rot.per}{percentage of vertical words}

  \item{max_min}{relative scales to adjust the size
  difference between largest and smallest word, by default
  the largest word is written with 10 times as large font
  than the smallest}

  \item{scale}{a fraction of the available space on figure
  that will be covered with the bounding boxes of words}

  \item{min.freq}{minimal frequency of words to be
  displayed}

  \item{max.words}{maximal number of words to be displayed}

  \item{random.order}{plot words in random order. If false,
  they will be plotted in decreasing frequency}

  \item{colors}{vector of colors fro the words. This vector
  will be extrapolated into as many colors as needed,
  starting with the first color for lower frequencies and
  ending with last color for higher frequencies.}

  \item{random.colors}{if true, assigns random color for
  the words.}

  \item{algorithm}{algorithm to find positions of words
  possible values: "circle", "leftside" and "rightside".}

  \item{tryfit}{if TRUE the algorithm checks if all words
  fit to the figure, if not it tries gradually smaller
  values of scale parameter until everything fits}

  \item{add}{if TRUE adds the picture to existing plot.}
}
\description{
  Uses the algorithm from wordcloud package to calculate
  the positions of the words. then uses grid graphics to
  plot the words on screen. The shape of the wordcloud
  depends on the shape of the plotting window
}
\examples{
plotWordcloud(c("Audi", "Volkswagen", "Opel", "Porsche", "Mercedez", "BMW"), 8:3)
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

