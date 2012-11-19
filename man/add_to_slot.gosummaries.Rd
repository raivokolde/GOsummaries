\name{is.gosummaries}
\alias{add_to_slot.gosummaries}
\alias{is.gosummaries}
\alias{print.gosummaries}
\title{Functions to interact with gosummaries object}
\usage{
  is.gosummaries(x)

  \method{print}{gosummaries} (x, ...)

  add_to_slot.gosummaries(x, slot, values)
}
\arguments{
  \item{gosummaries}{a gosummaries object}

  \item{\dots}{not used}

  \item{x}{a gosummaries object}

  \item{slot}{the component slot name to be filled (e.g
  Title, Percentage, etc.)}

  \item{values}{list of values where each element
  corresponds to one component}
}
\description{
  Functions to interact with gosummaries object
}
\examples{
data(gs_kmeans)

# Add new title to the components
gs_kmeans_new = add_to_slot.gosummaries(gs_kmeans, "Title", as.list(paste("K-means cluster", 1:length(gs_kmeans))))

print(gs_kmeans_new)
}
\author{
  Raivo Kolde <rkolde@gmail.com>
}

