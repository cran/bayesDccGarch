\name{logLikDccGarch}
\alias{logLikDccGarch}
\title{	The logarithm of likelihood function of DCC-GARCH(1,1) Model.}

\usage{
logLikDccGarch(mY, omega = rep(0.03, ncol(mY)), alpha = rep(0.03, ncol(mY)), 
	beta = rep(0.8, ncol(mY)), a = 0.03, b = 0.8, gamma = rep(1, ncol(mY)), 
	tail = 10, errorDist = 2)
}

\description{ 
 	Compute the logarithm of likelihood function of DCC-GARCH(1,1) Model if \code{mY} is a matrix or the logarithm of likelihood function of GARCH(1,1) Model if \code{mY} is numeric vector.
}

\arguments{
  \item{mY}{a matrix of the data (\eqn{n \times k}).}
  \item{omega}{a numeric vector (\eqn{k \times 1}) with the the values of \eqn{\omega_i} parameters. Default: \code{rep(0.03, ncol(mY))}.}
  \item{alpha}{a numeric vector (\eqn{k \times 1}) with the the values of \eqn{\alpha_i} parameters. Default: \code{rep(0.03, ncol(mY))}.}
  \item{beta}{a numeric vector (\eqn{k \times 1}) with the the values of \eqn{\beta_i} parameters. Default: \code{rep(0.80, ncol(mY))}.}
  \item{a}{a numeric value of the \eqn{a} parameter. Default: \code{0.03}.}
  \item{b}{a numeric value of the \eqn{b} parameter. Default: \code{0.8}.}
  \item{gamma}{a numeric vector (\eqn{k \times 1}) with the values of \eqn{\gamma_i} parameters. Default: \code{rep(1.0, ncol(mY))}.}
  \item{tail}{a numeric value of \eqn{\nu} parameter if \code{errorDist = 2} or of \eqn{\delta} parameter if \code{errorDist = 3}. If \code{errorDist = 1} so this arguments is no used.}
  \item{errorDist}{a probability distribution for errors. Use \code{errorDist=1} for \eqn{SSNorm}, \code{errorDist=2} for \eqn{SST} or \code{errorDist=3} for \eqn{SSGED}. Default: \code{2}.}
}

\details{
	The log-likelihood of the model GARCH(1,1) is computed if \code{mY} has just one column. The arguments \code{a} and \code{b} are not consider in this case.
}

\value{
Return a list with the elements:
\item{$H}{ a matrix where the lines are the \eqn{H_t} values for t=1,...,n.}
\item{$value}{the value of the logarithm of likelihood function.}
}

\references{Fioruci, J.A., Ehlers, R.S. Andrade Filho, M.G. \emph{Bayesian multivariate GARCH models with dynamic correlations and asymmetric error distributions}, Journal of Applied Statistics, 41(2), 320--331, 2014. \url{http://dx.doi.org/10.1080/02664763.2013.839635}.
}

\author{
	Jose Augusto Fioruci, Ricardo Sandes Ehlers and Francisco Louzada
}


\seealso{
	\code{\link{bayesDccGarch-package}}, \code{\link{bayesDccGarch}}
}

\examples{

data(DaxCacNik)

Dax = DaxCacNik[,1]

######  log-likelihood function of GARCH(1,1) model with SST innovations ####
logLikDccGarch(Dax, omega=0.03, alpha=0.03, beta=0.8, gamma=0.7)$value

######  log-likelihood function of DCC-GARCH(1,1) model with SST innovations ####
logLikDccGarch(DaxCacNik, beta=c(0.82,0.91,0.85), gamma=c(0.7, 1.3, 1.7), tail=10)$value

}

\keyword{ GARCH }
\keyword{ DCC-GARCH }
