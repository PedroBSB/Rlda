\name{plot}
\alias{plot.rlda}
\alias{plot,Rlda}
\title{plot method for rlda object}


\description{Plot a rlda object.
The \code{plot} function returns three plots based on Theta matrix, Phi matrix and log-likelihood. }

\usage{
  \method{plot}{rlda}(object,...)
}

\arguments{

  \item{object}{a \code{rlda} object created by \code{rlda.binomial}, \code{rlda.bernoulli} or \code{rlda.multinomial} function}
  \item{...}{ other arguments may be useful.}
}

\seealso{\code{\link{rlda.binomial}}, \code{\link{rlda.bernoulli}},\code{\link{rlda.multinomial}} }

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

\keyword{plot}   
\keyword{LDA}
\keyword{Fuzzy}
\keyword{Cluster}
