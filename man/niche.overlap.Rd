\name{niche.overlap}
\Rdversion{1.1}
\alias{niche.overlap}

\title{
Quantification of Niche Overlap
}

\description{
This function quantifies the degree of niche overlap using the statistics D and I (as proposed by Warren et al., 2008) based on Schoeners D (Schoener, 1968) and Hellinger Distances (van der Vaart, 1998). Niche overlap can be calculated either from the probability surfaces of entire ecological niche models (Warren et al., 2008) or from single dimensions of the climatic niches (Evans et al., 2009).
}

\usage{
niche.overlap(x)
}

\arguments{
  \item{x}{Might take one of the follwing forms: (1) a data frame containing predicted niche occupancy (PNO) profiles, e.g., as derived by \code{\link{pno}}; (2) a vector of filenames referencing to probability surfaces saved in ASCII-format with an ArcGIS-compatible header; (3) a list containing probability surfaces stored as objects of class \code{\link[sp]{SpatialGrid}}.
  }
}

\value{
A matrix of class \code{"niolap"}. The upper triangle contains pairwise comparisons of niche overlap in terms of D, whereas the lower triangle contains values of I.
}

\references{
Evans, M. E. K., S. A. Smith, R. S. Flynn, and M. J. Donoghue. 2009. Climate, niche evolution,
and diversification of the 'bird-cage evening primroses' (\emph{Oenothera}, sections \emph{Anogra} and 
\emph{Kleinia}). \emph{Am. Nat.} \bold{173}: 225-240.

Schoener, T.W. 1968. \emph{Anolis} lizards in Bimini: resource partitioning in a complex fauna. \emph{Ecology} \bold{49}: 704-726.

Van der Vaart, A.W. 1998. Asymptotic statistics. Cambridge University Press, U.K.

Warren, D. L., R. E. Glor, & M. Turelli. 2008. Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution} \bold{62}: 2868-2883.
}

\author{Christoph Heibl}

\seealso{\code{\link{pno}}, 
\code{\link{niche.equivalency.test}}, 
\code{\link{bg.similarity.test}},
\code{\link{age.range.correlation}}
}

\examples{
# load PNOs for Oxalis sect. Palmatifoliae
data(PNO)

# niche overlap on a annual precipitation gradient:
no <- niche.overlap(PNO$AnnualPrecipitation)

# upper triangle: based on Schoeners D
# lower triangle: based on Hellinger distances
print(no)	
}

