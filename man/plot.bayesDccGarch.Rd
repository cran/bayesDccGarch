\name{plot.bayesDccGarch}
\alias{plot.bayesDccGarch}
\alias{plot}
\title{Plotting volatilities for Bayesian DCC-GARCH model}

\usage{\method{plot}{bayesDccGarch}(x, ts.names=NULL, colors = c("grey","red"), \dots)
}

\arguments{
\item{x}{Object of class \dQuote{bayesDccGarch}.}
\item{ts.names}{a vector of length \eqn{k} with the names of the time series.}
\item{colors}{a vector with the colors for plotting the returns and volatilities.}
\item{\dots}{additional arguments for \code{plot} function}
}

\description{Produces a plot of time series and the volatilities. This is a particular case of \code{plotVol} function.}

\value{No return value}

\references{
Fioruci, J.A., Ehlers, R.S., Andrade Filho, M.G. \emph{Bayesian multivariate GARCH models with dynamic correlations and asymmetric error distributions}, Journal of Applied Statistics, 41(2), 320--331, 2014a. <doi:10.1080/02664763.2013.839635>

Fioruci, J.A., Ehlers, R.S., Louzada, F. \emph{BayesDccGarch - An Implementation of Multivariate GARCH DCC Models}, ArXiv e-prints, 2014b. \url{https://ui.adsabs.harvard.edu/abs/2014arXiv1412.2967F/abstract}.
}

\author{Ricardo Sandes Ehlers, Jose Augusto Fiorucci and Francisco Louzada}

\seealso{ \code{\link{bayesDccGarch-package}}, \code{\link{bayesDccGarch}}, \code{\link{plotVol}} }

\examples{

\donttest{
data(DaxCacNik)

mY = DaxCacNik

out = bayesDccGarch(mY, nSim=1000)
plot(out)
}

}

\keyword{plot}
\keyword{volatility}
