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
    \item{\code{aou530}}{Bird specie 530 presence (1) or absence (0)}
    \item{\code{aou540}}{Bird specie 540 presence (1) or absence (0)}
    \item{\code{aou590}}{Bird specie 590 presence (1) or absence (0)}
    \item ...
    \item{\code{aou4461}}{Bird specie 4461 presence (1) or absence (0)}
    \item{\code{aou5860}}{Bird specie 5860 presence (1) or absence (0)}
  }
}
\source{Pardieck, K.L., D.J. Ziolkowski Jr., M. Lutmerding, K. Campbell and M.-A.R. Hudson. 2017. North American Breeding Bird Survey Dataset 1966 - 2016, version 2016.0. U.S. Geological Survey, Patuxent Wildlife Research Center. \url{https://www.pwrc.usgs.gov/bbs/RawData/}doi:10.5066/F7W0944J.}
\examples{
data(birds)
}
\keyword{datasets}
