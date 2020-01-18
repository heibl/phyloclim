\name{phyloclim-package}
\alias{phyloclim-package}
\alias{phyloclim}
\docType{package}
\title{Phyloclimatic Data Analysis}
\description{
This package implements some recently developed methods in phyloclimatic modeling.
}
\details{
\tabular{ll}{
 Package: \tab phyloclim\cr
 Type: \tab Package\cr
Version: \tab 0.9.8\cr
Date: \tab 2020-01-18\cr
 License: \tab GPL (>= 2)\cr
}
\code{phyloclim} is intended to provide methods for the quantification of niche occupancy and niche overlap and combine this data with phylogenetic trees. 
Niche overlap can be quantified either from profiles of environmental suitability for a given environmental variable (PNO or predicted niche occupancy profiles, see \code{\link{pno}}) or from two-dimensional suitability surfaces (environmental niche models (ENMs) or species distribution models (SDMs)) with \code{\link{niche.overlap}}. Randomization tests for testing the significance of niche overlap in terms of niche identity and niche similarity are available in \code{\link{niche.equivalency.test}} and \code{\link{bg.similarity.test}}.
Phylogenetic data can be integrated into comparative methods such as \code{\link{age.range.correlation}} and the estimation of ancestral environmental tolerances (\code{\link{anc.clim}}), which is a non-parametric approach based on resampling of PNO profiles.
The functions in \code{phyloclim} are in an early stage of development, so they should be used with care! Any feedback is welcome!
}
\references{
Evans, M. E. K., S. A. Smith, R. S. Flynn, and M. J. Donoghue. 2009. Climate, niche evolution,
and diversification of the 'bird-cage evening primroses' (\emph{Oenothera}, sections \emph{Anogra} and 
\emph{Kleinia}). \emph{Am. Nat.} \bold{173}: 225-240.
Fitzpatrick, B.M & Turelli, M. 2006. The geography of mammalian speciation: mixed signals from phylogenies and range maps. \emph{Evolution} \bold{60}: 601-615.
Phillips, S.J, M. Dudik, & R.E. Schapire. 2006. Maximum entropy modeling of species geographic distributions. \emph{Ecological Modeling} \bold{190}: 231-259.
	
Warren, D., R.E. Glor, & M. Turelli. 2008. Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution} \bold{62}: 2868-2883.
}
\author{
Maintainer: Christoph Heibl (\email{christoph.heibl@gmx.net})
}
\keyword{package}
\seealso{
\code{\link[ape:ape-package]{ape}}
}
