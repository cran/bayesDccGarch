\name{bayesDccGarch-package}
\alias{bayesDccGarch-package}
\docType{package}

\title{
	bayesDccGARCH: Methods and tools for Bayesian analysis of DCC-GARCH(1,1) Model. 
}

\description{
	In this package we implemented functions for Bayesian analysis of DCC-GARCH(1,1) Model using the same modelling of Fioruci et al (2014). 
	Several probabilities distributions are available for the errors which can model both skewness and heavy tails.
}

\details{
	\tabular{ll}{
		Package: \tab bayesDccGarch\cr
		Type: \tab Package\cr
		Version: \tab 1.0\cr
		Date: \tab 2014-10-30\cr
		License: \tab GPL (>=2.14) \cr
    }
	
	\code{bayesDccGarch(mY, n_sim = 10000)}
}


\author{
	Jose Augusto Fioruci, Ricardo Sandes Ehlers and Francisco Louzada.
	Maintainer: Jose Augusto Fioruci <jafioruci@gmail.com> and Ricardo Sandes Ehlers <ehlers@icmc.usp.br>
}

\references{Fioruci, J.A., Ehlers, R.S. Andrade Filho, M.G. \emph{Bayesian multivariate GARCH models with dynamic correlations and asymmetric error distributions}, Journal of Applied Statistics, 41(2), 320--331, 2014. \url{http://dx.doi.org/10.1080/02664763.2013.839635}.
}

\keyword{ package }
\keyword{ multivariate }
\keyword{ DCC-GARCH }
\keyword{ GARCH }
\keyword{ Bayes }

\seealso{
Available functions:
\code{\link{bayesDccGarch}}, \code{\link{logLikDccGarch}}, \code{\link{dssnorm}}, \code{\link{dsst}}, \code{\link{dssged}}, \code{\link{plotVol}}

For classical estimation see:
\code{\link[ccgarch]{dcc.estimation}}
}

\examples{

data(DaxCacNik)

mY = DaxCacNik[1:10,] # more data is necessary
Dax = mY[,1]

### DCC-GARCH(1,1) ###
 out = bayesDccGarch(mY, nSim=1000)
 summary(out$MC)


### GARCH(1,1) ###
# out = bayesDccGarch(Dax)
# summary(out$MC) 
 
### DCC-GARCH(1,1) with SSGED innovations ###
# out = bayesDccGarch(mY, tail_ini=1.5, errorDist=3)
# summary(out$MC)
### Plotting volatilities #
# plotVol(mY, out$H[,c("H_1,1","H_2,2","H_3,3")], c("DAX","CAC40","NIKKEI"))


######  log-likelihood function of GARCH(1,1) model with SST innovations ####
logLikDccGarch(Dax, omega=0.03, alpha=0.03, beta=0.8, gamma=0.7)$value

######  log-likelihood function of DCC-GARCH(1,1) model with SST innovations ####
logLikDccGarch(mY, omega=c(0.02,0.01,0.03), alpha=c(0.02,0.01,0.03))$value

}
