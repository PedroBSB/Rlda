\name{summary}
\alias{summary}
\alias{summary.rlda}
\title{Provide summary information.}
\description{
Takes a \code{rlda} object produced by \code{rlda.binomial}, \code{rlda.bernoulli} or \code{rlda.multinomial} and produces a summary from it.
}
\usage{
  \method{summary}{rlda}(object,...)
}
\arguments{
  \item{object}{a \code{rlda} object as produced by \code{rlda.binomial}, \code{rlda.bernoulli} or \code{rlda.multinomial}.}
  \item{...}{ other arguments may be useful.}
}

\details{
Summarizes the Gibbs Samping results and arguments.
}

\author{
  \itemize{
    \item
      Pedro Albuquerque.\cr
      \email{pedroa@unb.br}\cr
      \url{http://pedrounb.blogspot.com/}
   \item
      Denis Valle.\cr
      \email{drvalle@ufl.edu}\cr
      \url{http://denisvalle.weebly.com/}
   \item
      Daijiang Li.\cr
      \email{daijianglee@gmail.com}\cr
  }
}

\keyword{summary}
\keyword{LDA}
\keyword{Fuzzy}
\keyword{Cluster}

\seealso{\code{\link{rlda.binomial}}, \code{\link{rlda.bernoulli}},\code{\link{rlda.multinomial}} }

\examples{
\dontrun{
if(requireNamespace("INLA", quietly = TRUE)){
  require("INLA", quietly = TRUE)
  data(Catheter)
  res <- rlda(data = Catheter)
  summary(res)
  }
}
}