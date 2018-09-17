\name{pno.weighted.mean}
\Rdversion{1.1}
\alias{pno.weighted.mean}

\title{
Weighted means of niche dimensions
}

\description{
This function calculates the mean of an environmental variable weighted by the cumulative probability of each value of the environmental gradient extracted from a species distibution model (e.g. MAXENT)
}

\usage{
pno.weighted.mean(x, subset = NULL, normalize = TRUE)
}

\arguments{
  \item{x}{A data frame or matrix with columns corresponding to species and rows corresponding to values along an environmental gradient. The first columns contains the environmental variable, the remaining colums probabilities of suitability.
  }
  \item{subset}{A vector of mode \code{"character"} which can be used to restrict the calculation of weighted means to those columsn in \code{x} whose column names match \code{subset}; defaults to \code{NULL}. }
  
  \item{normalize}{A logical, indicating whether to normalize PNOs prior to calculation of weighted means.}
}

\value{
A vector of weighted means corresponding to the columns in \code{x}.
}

\references{
Evans, M. E. K., S. A. Smith, R. S. Flynn, and M. J. Donoghue. 2009. Climate, niche evolution,
and diversification of the 'bird-cage evening primroses' (\emph{Oenothera}, sections \emph{Anogra} and 
\emph{Kleinia}) \emph{Am. Nat.} \bold{173}: 225-240.
}

\author{Christoph Heibl}

\seealso{
\code{\link{pno}}, 
\code{\link{anc.clim}}, 
\code{\link{plotAncClim}}
}