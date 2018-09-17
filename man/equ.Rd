\name{enneaphylla_laciniata}
\Rdversion{1.1}
\alias{equ}
\alias{sim}
\docType{data}
\title{
Niche overlap for Oxalis enneaphylla and O. lacinita
}
\description{
Objects \code{equ} and \code{sim} contain the results of the \code{\link{niche.equivalency.test}} and the \code{\link{bg.similarity.test}} (Warren et al. 2008) for \emph{Oxalis enneaphylla} and \emph{O. laciniata} from section \emph{Palmatifoliae}.
}

\usage{
data(equ)
data(sim)
}

\format{
A list containing five elements:

\describe{
 \item{method}{Name of the test.
 }
 \item{species}{Names of the two species compared.
 }
 \item{statistic}{Test statistics D and I of niche overlap based on Schoeners D and modified Hellinger distances.
 }
 \item{p.value}{p-values for the test statistics (if \code{method = "niche equivalency test"}).
 }
 \item{null}{Formulation of the null hypothesis.
 }
 \item{null.distribution}{Null distributions of D and I derived from randomization (if \code{method = "niche equivalency test"}).
 }
\item{ci.x.randomY}{Confidence intervals based on null distributions (if \code{method = "background similarity test"}).
}
\item{ci.y.randomX}{Confidence intervals based on null distributions (if \code{method = "background similarity test"}).
}
\item{nd.x.randomY}{Null distributions of D and I derived from randomization (if \code{method = "background similarity test"}).
}
\item{nd.y.randomX}{Null distributions of D and I derived from randomization (if \code{method = "background similarity test"}).
}
}
}

\source{
C. Heibl, unpublished data.
}

\references{
Warren, D., R.E. Glor, & M. Turelli. 2008. Environmental niche equivalency versus conservatism: quantitative approaches to niche evolution. \emph{Evolution}. \bold{62}: 2868-2883.
}

\examples{
# plot results of niche equivalency test
data(equ)
plot(equ)

# plot results of background similarity test
data(sim)
plot(sim)
}
\keyword{datasets}
