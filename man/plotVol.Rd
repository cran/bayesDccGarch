\name{plotVol}
\alias{plotVol}

\title{
  Plotting volatilities of time series
}

\description{
 Plotting method for volatilities of time series.
}

\usage{
plotVol(mY, vol, ts.names=paste("TS_", 1:ncol(mY), sep=""), colors = c("grey","red"), \dots)
}

\arguments{
  \item{mY}{a matrix of the data (\eqn{n \times k}).}
  \item{vol}{a matrix (\eqn{n \times k}) with the volatility estimates.}
  \item{ts.names}{a vector of length \eqn{k} with the names of the time series.}
  \item{colors}{a vector with name of the colors for plotting the returns and volatilities.}
  \item{\dots}{additional arguments for \code{plot} function}
}

\references{
Fioruci, J.A., Ehlers, R.S., Andrade Filho, M.G. \emph{Bayesian multivariate GARCH models with dynamic correlations and asymmetric error distributions}, Journal of Applied Statistics, 41(2), 320--331, 2014a. \url{http://dx.doi.org/10.1080/02664763.2013.839635}.

Fioruci, J.A., Ehlers, R.S., Louzada, F. \emph{BayesDccGarch - An Implementation of Multivariate GARCH DCC Models}, ArXiv e-prints, 2014b. \url{http://adsabs.harvard.edu/abs/2014arXiv1412.2967F}.
}

\author{Ricardo Sandes Ehlers, Jose Augusto Fioruci and Francisco Louzada}

\seealso{ \code{\link{bayesDccGarch-package}}, \code{\link{bayesDccGarch}} }

\examples{

data(DaxCacNik)

mY = DaxCacNik[1:10,] # more data is necessary

out = bayesDccGarch(mY, nSim=1000)
plotVol(mY, out$H[,c("H_1,1","H_2,2","H_3,3")], c("DAX","CAC40","NIKKEI"))

}

\keyword{plot}
\keyword{volatility}


