\name{birds}
\alias{birds}
\docType{data}
\title{
Breeding Bird Survey
}
\description{
This dataset is a subset of the Breeding Bird Survey.
}
\usage{data("birds")}
\format{
  A data frame with 13608 observations and 384 variables.
  \describe{
    \item{\code{loc.id}}{Location index}
    \item{\code{aouXXXX}}{Bird specie presence (1) or absence (0)}
  }
}
\source{\url{https://www.pwrc.usgs.gov/bbs/RawData/}}
\examples{
data(birds)
}
\keyword{datasets}
