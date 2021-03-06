% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bg.similarity.test.R,
%   R/niche.equivalency.test.R, R/niche.tests.R, R/plot.ntest.R, R/print.ntest.R
\name{bg.similarity.test}
\alias{bg.similarity.test}
\alias{niche.equivalency.test}
\alias{niche.tests}
\alias{plot.ntest}
\alias{print.ntest}
\title{Niche Equivalency and Background Similarity Test}
\usage{
bg.similarity.test(
  p,
  env,
  n = 99,
  study.area.y = "mcp",
  conf.level = 0.95,
  app,
  dir
)

niche.equivalency.test(p, env, n = 99, app, dir)

\method{plot}{ntest}(x, ...)

\method{print}{ntest}(x, ...)
}
\arguments{
\item{p}{A \code{\link[sp]{SpatialPointsDataFrame}} or a simple data frame
containing the \bold{p}resence points. In the latter case the first column
contains the species names, the second and third column longitude and
latitude (see SWD-formatted (=Samples-With-Data) files in the MAXENT
tutorial).}

\item{env}{An object of class \code{\link[sp]{SpatialGridDataFrame}}
containing the environmental covariates.}

\item{n}{An integer giving the number of permutations of the original data
(default: \code{n = 99}).}

\item{study.area.y}{Defines the study area of the second species Y (see
Warren et al. 2008). Can be \code{"env"} (study area corresponds to
\code{env}) or \code{"mcp"} (study area is the area of a convex hull around
the presence points of the second species Y).}

\item{conf.level}{A real number between 0 and 1 setting the confidence level
of the confidence intervals to be calculated.}

\item{app}{A character string giving the path to the MAXENT application.}

\item{dir}{A character string giving the name of a directory where the input
and output data for MAXENT will be saved. Already existing directories will
be overwritten \bold{without} a warning. If \code{dir} is left empty the
data will be written to a temporary directory, which will be deleted after
execution.}

\item{x}{An object of class \code{ntest}.}

\item{\dots}{Further arguments passed to or from other methods.}
}
\value{
\code{niche.equivalency.test} gives a list with six elements:
  \item{method}{Name of the test} 
  \item{species}{Names of the two species
  compared} 
  \item{null}{Formulation of the null hypothesis}
  \item{statistic}{Statistics of niche overlap D based on Schoeners D and
  modified Hellinger distances} 
  \item{p.value}{p-values associated with the
  statistics} 
  \item{null.distribution}{Null distributions of D and I derived
  from randomization}

\code{bg.similarity.test} gives a list with eight elements:
  \item{method}{Name of the test} \item{species}{Names of the two species
  compared} \item{null}{Formulation of the null hypothesis}
  \item{statistic}{Statistics of niche overlap D based on Schoeners D and
  modified Hellinger distances} 
  \item{ci.x.randomY}{Confidence interval for D and I based on the comparison
  of the first species against a set of random presence points from the study
  area of the second species}
  \item{ci.y.randomX}{Confidence interval for D and I based on the comparison
  of the second species against a set of random presence points from the
  study area of the first species}
  \item{nd.x.randomY}{Null distributions of D and I calculated from the
  comparison of the first species against a set of random presence points
  from the study area of the second species}
  \item{nd.y.randomX}{Null distributions of D and I calculated from the
  comparison of the second species against a set of random presence points
  from the study area of the first species}
}
\description{
Hypothesis testing as proposed by
  \insertCite{warrenturelli2008;textual}{phyloclim} based on the generation of
  pseudoreplicate datasets. The \strong{niche equivalency (or identity) test} asks
  whether the ecological niche models (ENMs) of two species are more
  different than expected if they were drawn from the same underlying
  distribution. The \strong{background similarity test} asks whether ENMs drawn from
  populations with partially or entirely non-overlapping distributions are
  any more different from one another than expected by chance.
}
\details{
An installation of MAXENT \insertCite{phillipsdudik2008}{phyloclim}
  is required in order to run \code{niche.equivalency.test} and
  \code{bg.similarity.test}. Both functions use the logistic output of MAXENT
  estimated using auto features.

  By default, the environmental covariates given with \code{env} are assumend
  to be \emph{continuous}. In order to use \emph{categorical} environmental
  covariates, you have to prepend \code{"cat_"} to the layer name, e.g.
  \code{"cat_landuse"}.
}
\note{
These functions have been tested with MAXENT 3.3.4
}
\examples{
# path to MAXENT
# --------------
maxent.exe <- file.path(system.file(package="dismo"), "java/maxent.jar")

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
\dontrun{
  if (file.exists(maxent.exe)){
    net <- niche.equivalency.test(samples, preds, reps, maxent.exe)
    net; plot(net)
    bst <- bg.similarity.test(samples, preds, reps, app = maxent.exe)
    bst; plot(bst)
  } else {
    message("get a copy of MAXENT (see Details)")
  }
}
}
\references{
\insertRef{phillipsdudik2008}{phyloclim}

MAXENT webseite: \url{https://biodiversityinformatics.amnh.org/open_source/maxent/}

\insertRef{warrenturelli2008}{phyloclim}

\insertRef{warrenseifert2011}{phyloclim}
}
