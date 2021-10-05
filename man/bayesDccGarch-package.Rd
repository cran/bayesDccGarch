\name{bayesDccGarch-package}
\alias{bayesDccGarch-package}
\docType{package}

\title{
	bayesDccGARCH: Methods and tools for Bayesian analysis of DCC-GARCH(1,1) Model.
}

\description{
	In this package we implemented functions for Bayesian analysis of DCC-GARCH(1,1) Model using the same modelling of Fioruci et al (2014a).
	Several probabilities distributions are available for the errors which can model both skewness and heavy tails.
	See Fioruci et al (2014b) for more details about the package.
}

\details{
	\tabular{ll}{
		Package: \tab bayesDccGarch\cr
		Type: \tab Package\cr
		Version: \tab 3.0.3\cr
		Date: \tab 2021-10-02\cr
		License: \tab GPL (>=2.14) \cr
    }

	\code{bayesDccGarch(mY, n_sim = 10000)}
}


\author{
	Jose Augusto Fiorucci, Ricardo Sandes Ehlers and Francisco Louzada.
	Maintainer: Jose Augusto Fiorucci <jafiorucci@gmail.com>
}

\references{
Fioruci, J.A., Ehlers, R.S., Andrade Filho, M.G. \emph{Bayesian multivariate GARCH models with dynamic correlations and asymmetric error distributions}, Journal of Applied Statistics, 41(2), 320--331, 2014a, <doi:10.1080/02664763.2013.839635>.

Fioruci, J.A., Ehlers, R.S., Louzada, F. \emph{BayesDccGarch - An Implementation of Multivariate GARCH DCC Models}, ArXiv e-prints, 2014b. \url{https://ui.adsabs.harvard.edu/abs/2014arXiv1412.2967F/abstract}
}

\keyword{ package }
\keyword{ multivariate }
\keyword{ DCC-GARCH }
\keyword{ GARCH }
\keyword{ Bayes }

\seealso{
Available functions:
\code{\link{bayesDccGarch}}, \code{\link{update}}, \code{\link{predict}}, \code{\link{plot}}, \code{\link{logLikDccGarch}}, \code{\link{dssnorm}}, \code{\link{dsst}}, \code{\link{dssged}}, \code{\link{plotVol}}
}

\examples{

\donttest{
data(DaxCacNik)

out = bayesDccGarch(DaxCacNik)

summary(out)

plot(out)

}

}
