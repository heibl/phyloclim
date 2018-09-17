\name{plotPNO}
\Rdversion{1.1}
\alias{plotPNO}

\title{
Plot Predicted Niche Occupancy Profiles
}

\description{
This function can be used to plot predicted niche occupancy profiles (PNOs). PNOs can be obtained in a geographical information system by summing the cumulative probabilies of each climatical value for a species distribution model (SDM).
}

\usage{
plotPNO(x, subset = NULL, thinning = NULL, xlab = NULL, 
	tail_threshold = 0, wm = FALSE, legend.pos = "topleft")
}

\arguments{
  \item{x}{A data frame or matrix with columns corresponding to species and rows corresponding to values along an environmental gradient. The first columns contains the environmental variable, the remaining colums probabilities of suitability.
}
  \item{subset}{A vector of mode \code{"character"} which can be used to restrict the calculation of weighted means to those columsn in \code{x} whose column names match \code{subset}; defaults to \code{NULL}.}
  
  \item{thinning}{An integer that can be used to thin fuzzy PNOs prior to plotting; defaults to \code{NULL}.}
  
  \item{xlab}{A character string given the label for the x-axis.}
  
  \item{tail_threshold}{A numeric that can be used cut long tails of PNOs; defaults to \code{0}.}
  
  \item{wm}{A logical indicating if weighted mean should added for each species.}
  \item{legend.pos}{Controls the position of the legend. Might eihter be a list object containing x and y coordinates (such as e.g. returned by \code{\link{locator}}) of the \bold{topleft corner} of the legend box or one of the following: \code{"topleft"} (default), \code{"bottomleft"}, \code{"topright"}, or \code{"bottomright"}. If \code{legend.pos == NULL} the plotting of the legend is suppressed.}
}

\references{
Evans, M. E. K., S. A. Smith, R. S. Flynn, and M. J. Donoghue. 2009. Climate, niche evolution,
and diversification of the 'bird-cage evening primroses' (\emph{Oenothera}, sections \emph{Anogra} and 
\emph{Kleinia}). \emph{Am. Nat.} \bold{173}: 225-240.
}

\author{Christoph Heibl}

\seealso{
\code{\link{pno}}
}
\examples{
# load PNOs for Oxalis sect. Palmatifoliae
data(PNO)

# plot predicted niche occupany for annual mean temperature
plotPNO(x = PNO$AnnualMeanTemperature, 
	xlab = "Annual Mean Temperature (degree C)")
	
# same plot, but with weighted means added
plotPNO(x = PNO$AnnualMeanTemperature, 
	xlab = "Annual Mean Temperature (degree C)", wm = TRUE)
}
