\name{hypothesis.testing}
\Rdversion{1.1}
\alias{niche.equivalency.test}
\alias{niche.identity.test}
\alias{identity.test}
\alias{bg.similarity.test}
\alias{print.ntest}
\alias{plot.ntest}

\title{Niche Equivalency and Background Similarity Test}

\description{
Hypothesis testing as proposed by Warren et al. (2008) based on the generation of pseudoreplicate datasets. The niche equivalency (or identity) test asks whether the ecological niche models (ENMs) of two species are more different than expected if they were drawn from the same underlying distribution. The background similarity test asks whether ENMs drawn from populations with partially or entirely non-overlapping distributions are any more different from one another than expected by random chance. 
}

\usage{
niche.equivalency.test(p, env, n = 99, app, dir)

bg.similarity.test(p, env, n = 99, conf.level = .95, app, dir)


\method{print}{ntest}(x, \dots)

\method{plot}{ntest}(x, \dots)
}

\arguments{
  \item{p}{a \code{\link[sp]{SpatialPointsDataFrame}} or a simple data frame containing the \bold{p}resence points. In the latter case the first column contains the species names, the second and third column longitude and latitude (see SWD-formatted (=Samples-With-Data) files in the MAXENT tutorial).}
  
  \item{env}{an object of class \code{\link[sp]{SpatialGridDataFrame}} containing the environmental covariates.}
  
  \item{n}{an integer giving the number of permutations of the original data (default: \code{n = 99}).}
  
  \item{conf.level}{a real number between 0 and 1 setting the confidence level of the confidence intervals to be calculated.}
  
  \item{app}{a character string giving the path to the MAXENT application.}
  
  \item{dir}{a character string giving the name of a directory where the input and output data for MAXENT will be saved. Already existing directories will be overwritten \bold{without} a warning. If \code{dir} is left empty the data will be written to a temporary directory, which will be deleted after execution.}
  
  \item{x}{an object of class \code{ntest}.}

  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
An installation of MAXENT (Phillips et al., 2006; \url{http://www.cs.princeton.edu/~schapire/maxent/}) is required in order to run \code{niche.equivalency.test} and \code{bg.similarity.test}. Both functions use the logistic output of MAXENT estimated using auto features. 

By default, the environmental covariates given with \code{env} are assumend to be \emph{continuous}. In order to use \emph{categorical} environmental covariates, you have to prepend \code{"cat_"} to the layer name, e.g. \code{"cat_landuse"}.
}

\value{
\code{niche.equivalency.test} gives a list with six elements:
\item{method}{name of the test.}
\item{species}{names of the two species compared.}
\item{null}{formulation of the null hypothesis.}
\item{statistic}{statistics of niche overlap D based on Schoeners D and modified Hellinger distances.}
\item{p.value}{p-values associated with the statistics.}
\item{null.distribution}{null distributions of D and I derived from randomization.}

\code{bg.similarity.test} gives a list with eight elements:
\item{method}{name of the test.}
\item{species}{names of the two species compared.}
\item{null}{formulation of the null hypothesis.}
\item{statistic}{statistics of niche overlap D based on Schoeners D and modified Hellinger distances.}
\item{ci.x.randomY}{confidence interval for D and I based on the comparison of the first species against a randomized background derived from the second species.}
    \item{ci.y.randomX}{confidence interval for D and I based on the comparison of the second species against a randomized background derived from the first species.}
    \item{nd.x.randomY}{null distributions of D and I calculated from the comparison of the first species against a randomized background derived from the second species.}
\item{nd.y.randomX}{null distributions of D and I calculated from the comparison of the second species against a randomized background derived from the first species.}

}

\references{	
Phillips, S.J, M. Dudik, & R.E. Schapire. 2006. Maximum entropy modeling of species geographic distributions. \emph{Ecological Modeling} \bold{190}: 231-259.	
	
Warren, D., R.E. Glor, & M. Turelli. 2008. Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution}. \bold{62}: 2868-2883.
}

\author{Christoph Heibl}

\note{These functions have been completely rewritten and have been tested with MAXENT 3.3.3k}

\seealso{\code{\link{niche.overlap}}}
\examples{
# path to MAXENT
# --------------
maxent.exe <- paste(system.file(package="dismo"), 
"/java/maxent.jar", sep = "")

# a data frame of coordinates where two species 
# have been detected ('presence points') and
# a raster stack of environmental covariables
# --------------------------------------
species <- c("enneaphylla", "laciniata")
data(sites)
samples <- sites[grep(paste(species, collapse = "|"), sites$spec), ]
data.path <- system.file("extdata", package = "phyloclim")
preds <- list.files(path = data.path, pattern = "[.]asc")
preds <- paste(data.path, preds, sep = "/")
preds <- stack(lapply(X = preds, FUN = raster))

# testing against 9 permutations of the data
# -------------------------------------------
reps <- 9

# run hypothesis tests
# --------------------
if (file.exists(maxent.exe)){
  net <- niche.equivalency.test(samples, preds, reps, maxent.exe)
  net; plot(net)
  bst <- bg.similarity.test(samples, preds, reps, app = maxent.exe)
  bst; plot(bst)
} else {
  message("get a copy of MAXENT (see Details)")
}
}