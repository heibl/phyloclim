\name{palmatifoliae}
\Rdversion{1.1}
\alias{palmatifoliae}
\alias{PNO}
\alias{sites}
\alias{tree}
\docType{data}
\title{
\emph{Oxalis} section \emph{Palmatifoliae}
}
\description{
	This is a data set of \emph{Oxalis} sect. \emph{Palmatifoliae} derived from a larger scale study by Heibl & Renner (2012). The \code{sites} data frame gives the coordinates (WGS1984) of the sites where \emph{Oxalis} sect. Palmatifoliae has been collected. These were taken as presence points in distribution modeling and the list of predicted niche occupancy profiles (\code{PNO}). Finally, \code{tree} holds a phylogenetic hypothesis for the section extracted from a relaxed molecular clock model of the phylogeny of Oxalidales based approx. 7500 bp of chloroplast and nuclear markers.
For 19 bioclimatic variables summarizing temperature and precipiation and five species of \emph{Oxalis}, the suitability or probability of occurence is given along the environmental gradients. 
}
\usage{
	data(PNO)
	data(sites)
	data(tree)
	}
\format{
  \describe{
    \item{\code{sites}}{a data frame representing 290 collections of \emph{Oxalis} sect. \emph{Palmatifoliae}.}
    \item{\code{PNO}}{a list consisting of 19 matrixes. Each matrix corresponds to an environmental variable and contains the predicted niche occupancy profile (PNO) for each species of \emph{Oxalis} sect. \emph{Palmatifoliae}, whereby column 1 of each matrix gives the values of the environmental variable and columns 2-6 give probabilites of occurance of each species along the environmental gradient in column 1.}
    \item{\code{tree}}{an object of class \code{"phylo"} (see Value section of \code{\link[ape:read.tree]{read.tree}}).}
  }
}
\source{
Heibl, C. & S.S. Renner. Arid-adapted Distribution models and a dated phylogeny for Chilean \emph{Oxalis} species reveal occupation of new habitats by different lineages, not rapid adaptive radiation. 2012. \emph{Syst. Biol.} \bold{61}(5): 823-834.
}
\examples{
# 290 collections of Palmatifoliae
# --------------------------------
data(sites)
plot(sites$long, sites$lat,
  xlab = "Longitude", ylab = "Latitude")

## PNO profiles along 19 bioclimatic variables
## -------------------------------------------
data(PNO); names(PNO)
temp <- names(PNO)[1]
plotPNO(PNO[[temp]], xlab = temp)

## phylogenetic hypothesis for Palmatifoliae
## -----------------------------------------
data(tree)
plot(tree)
}
\keyword{datasets}
