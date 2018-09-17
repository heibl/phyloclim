\name{age.range.correlation}
\Rdversion{1.1}
\alias{age.range.correlation}

\title{Age-Range Correlation
}

\description{This function can be used to test for phylogenetic signal in patterns of niche overlap (Warren et al., 2008) based on the age-range correlation (ARC) as implemented by Turelli & Fitzpatrick (2006).  
}

\usage{
age.range.correlation(phy, overlap, tri = "upper", n = 10000)
}

\arguments{
  \item{phy}{an object of class \code{"phylo"}, which must be ultrametric.
}
  \item{overlap}{a matrix of class \code{"niolap"}, containing pairwise measures of niche overlap as produced by \code{\link{niche.overlap}}.
}

  \item{tri}{either \code{"upper"} or \code{"lower"}, indicating which triangle of the overlap matrix to use (see \code{\link{niche.overlap}}).
}
  \item{n}{an integer giving the number of iterations for the Monte Carlo resampling procedure.
}
}
\details{
In age-range correlation (ARC; \emph{sensu} Fitzpatrick & Turelli 2006), phylogentic relationships are taken into account when computing average range (or niche) overlap between two clades at node i by 

\deqn{
	\bar{o}_{i} = \sum_{j \in Clade1} \sum_{k \in Clade2} (\frac{1}{2})^{n_{jk} - 1} o_{jk}
	}{
	o(i) = sum( over j) sum(over k)  0.5^(n(ik) - 1) * o(jk)
	} 
	
where the double sum is over all species in the two clades, o(jk) denotes the overlap between species j and k, and n(jk) is the number of nodes separating the two species on the phylogeny. Plotting \eqn{\bar{o}_{i}}{o(i)} against node ages depicts the ARC for a given clade. Slopes and intercept derived from a linear model can be used to characterize speciation mode (allopatric versus sympatric) or niche evolution (conservatism versus flexibility) in the clade (see Fitzpatrick & Turelli [2006] and Warren et al. [2008], respectively).

The significance of the ARC is assessed via Monte Carlo simulation by random permutation of the overlap matrix to estimate the distribution  of slope and intercept unter the null hypothesis of no phylogenetic signal.
}

\value{
A list of four elements: 
\item{age.range.correlation}{A matrix with age and average overlap for each node.}
\item{linear.regression}{A linear model of mean niche (or range)overlap versus node age.}
\item{sig}{Significance of intercept and slope derived by randomization. \code{f} is the fraction of randomized datasets with intercept and slopes greater than the observed data. \code{p}-values are calculated as 2 * min(f, 1 - f).}
\item{MonteCarlo.replicates}{A matrix containing intercepts and slopes of randomized datasets.}
}

\references{
Fitzpatrick, B.M & Turelli, M. 2006. The geography of mammalian speciation: mixed signals from phylogenies and range maps. \emph{Evolution} \bold{60}: 601-615.	
	
Warren, D., R.E. Glor, & M. Turelli. 2008. Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution} \bold{62}: 2868-2883.
}

\author{Christoph Heibl}

\seealso{\code{\link{niche.overlap}}}

\examples{
# load PNOs for Oxalis sect. Palmatifoliae ...
data(PNO)

# ... and calculate niche overlap between species
no <- niche.overlap(PNO$AnnualMeanTemperature)

# load phylogeny and PNOs of Oxalis sect. Palmatifoliae
data(tree)

# age-range correlation
x <- age.range.correlation(phy = tree, overlap = no, n = 100)

# plot average niche overlap versus node age
plot(x$age.range.correlation)

# add a regression line
abline(x$linear.regression$coefficients)

# add regression lines from Monte Carlo randomization
apply(x$MonteCarlo.replicates, 1, abline, lwd = 0.2, col = "grey50")
}

