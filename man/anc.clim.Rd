\name{anc.clim}
\Rdversion{1.1}
\alias{anc.clim}

\title{
Estimation of Ancestral Climatic Tolerances
}

\description{
This function implements a nonparametric approach to estimate ancestral climatic tolerances proposed by Evans et al. (2009), using \code{\link[ape:ace]{ace}} from the \code{ape} package. 
}

\usage{
anc.clim(target, posterior = NULL, pno, n = 100, method = "GLS")
}

\arguments{
  \item{target}{An object of class \code{"phylo"}; a best estimate  (e.g., ML estimate) of topology on which to summarize mean estimated climatic tolerances, must be ultrametric.
}

  \item{posterior}{An object of class \code{"multiPhylo"} containing a sample of ultrametric trees from the posterior distribution of a Bayesian tree search to be used for reconstruction of ancestral climatic tolerances. \code{posterior = NULL} (the default) means that reconstruction uses \code{target} instead.
}
  \item{pno}{A data frame containing predicted niche occupancy (PNO) profiles, produced e.g. with \code{\link{pno}}.}
  
  \item{n}{An integer: The PNO profile given by \code{pno} argument will be resamples \code{n} times.}
  
   \item{method}{A character string specifying the method used for estimation. Two choices are possible: \code{"ML"} or \code{"GLS"} (see \code{\link[ape:ace]{ace}} for details).
}
}

\value{
A list containing three elements:
\item{tree}{An object of class \code{"phylo"}, the \code{target} tree used in the estimation.}
\item{means}{A vector containing the means of estimated ancestral climatic tolerances. The elements of \code{means} are in the same order as node numbers in element \code{tree}.}
\item{central.density}{A matrix, which for each of the terminal nodes in \code{tree} contains the minimum and the maximum of the 80 percent density interval of climatic tolerances (as sampled from the PNO).}

}

\references{
Evans, M. E. K., S. A. Smith, R. S. Flynn, and M. J. Donoghue. 2009. Climate, niche evolution,
and diversification of the 'bird-cage evening primroses' (\emph{Oenothera}, sections \emph{Anogra} and 
\emph{Kleinia}). \emph{Am. Nat.} \bold{173}: 225-240.
}

\author{Christoph Heibl}

\seealso{
\code{\link{pno}}, 
\code{\link{plotAncClim}}
}

\examples{
# load phylogeny and PNOs of Oxalis sect. Palmatifoliae
data(tree)
data(PNO)

# choose summer precipitation for analysis
clim <- PNO$PrecipitationWarmestQuarter

# estimate ancestral tolerances
ac <- anc.clim(target = tree, pno = clim, n = 100)

# visualize results
plotAncClim(ac, ylab = "Precipitation of warmest quarter (mm)")
}

