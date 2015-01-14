\name{densityFunctions}
\alias{dssnorm}
\alias{dsst}
\alias{dssged}

\title{Density functions of multivariate Standard Skew Norm, t-Student and GED distributions}

\description{ 
 	Compute the density function of Standard Skew Normal distribution (SSNORM) or 
	density function of Standard Skew t-Student distribution (SST) or 
	density function of Standard Skew GED distribution (SSGED)
	}

\usage{
  dssnorm(x, gamma=rep(1,length(x)), log=FALSE)
  
  dsst(x, gamma=rep(1,length(x)), nu=10, log=FALSE)
  
  dssged(x, gamma=rep(1,length(x)), delta=2, log=FALSE)
}

\arguments{
  \item{x}{a numeric vector for the point which the density will be computed.}
  \item{gamma}{a numeric vector for skew parameters. Must be positive.}
  \item{nu}{a numeric value of shape parameter of the multivariate Standard Skew t-Student distribution. Must be greater than 2.}
  \item{delta}{a numeric value of shape parameter of GED distribution. Must be positive.}
  \item{log}{logical; if TRUE, densities \code{p} are returned as \code{log(p)}.}
}

\value{Returns the computed value of the density.}

\references{Fioruci, J.A., Ehlers, R.S. Andrade Filho, M.G. \emph{Bayesian multivariate GARCH models with dynamic correlations and asymmetric error distributions}, Journal of Applied Statistics, 41(2), 320--331, 2014. \url{http://dx.doi.org/10.1080/02664763.2013.839635}.
}

\author{Jose Augusto Fioruci, Ricardo Sandes Ehlers and Francisco Louzada}

\seealso{ \code{\link{bayesDccGarch-package}} }


\examples{
### Univariate symmetric standard norm distributions ###
dssnorm(x=0)
dsst(x=0, nu=100)
dssged(x=0, delta=2)

### Univariate standard skew norm distributions ###
dssnorm(x=0, gamma=1.5)
dsst(x=0, gamma=1.5, nu=100)
dssged(x=0, gamma=1.5, delta=2)

### Multivariate standard skew norm distributions ###
dssnorm(x=c(0,0), gamma=c(1.5,0.7))
dsst(x=c(0,0), gamma=c(1.5,0.7), nu=100)
dssged(x=c(0,0), gamma=c(1.5,0.7), delta=2)
}

\keyword{ distribution }
\keyword{ density-function }
\keyword{ multivariate-distribution }
\keyword{ skew-distribution }