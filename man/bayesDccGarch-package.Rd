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
		Version: \tab 2.0\cr
		Date: \tab 2016-01-29\cr
		License: \tab GPL (>=2.14) \cr
    }
	
	\code{bayesDccGarch(mY, n_sim = 10000)}
}


\author{
	Jose Augusto Fiorucci, Ricardo Sandes Ehlers and Francisco Louzada.
	Maintainer: Jose Augusto Fiorucci <jafioruci@gmail.com> and Ricardo Sandes Ehlers <ehlers@icmc.usp.br>
}

\references{
Fioruci, J.A., Ehlers, R.S., Andrade Filho, M.G. \emph{Bayesian multivariate GARCH models with dynamic correlations and asymmetric error distributions}, Journal of Applied Statistics, 41(2), 320--331, 2014a. \url{http://dx.doi.org/10.1080/02664763.2013.839635}.

Fioruci, J.A., Ehlers, R.S., Louzada, F. \emph{BayesDccGarch - An Implementation of Multivariate GARCH DCC Models}, ArXiv e-prints, 2014b. \url{http://adsabs.harvard.edu/abs/2014arXiv1412.2967F}.
}

\keyword{ package }
\keyword{ multivariate }
\keyword{ DCC-GARCH }
\keyword{ GARCH }
\keyword{ Bayes }

\seealso{
Available functions:
\code{\link{bayesDccGarch}}, \code{\link{logLikDccGarch}}, \code{\link{dssnorm}}, \code{\link{dsst}}, \code{\link{dssged}}, \code{\link{plot}}, \code{\link{plotVol}}

For classical estimation see \emph{Package ccgarch}.
}

\examples{

data(DaxCacNik)  
Dax = DaxCacNik[,1]

### DCC-GARCH(1,1) ###
out = bayesDccGarch(DaxCacNik[1:10,], nSim=1000) # more data is necessary
summary(out)
plot(out) 

######## using more data and more simulation ###################
### GARCH(1,1) ###
# out = bayesDccGarch(Dax, nSim=20000)
# out2 = window(out, start=5000) # removing the first 5000 as burning
# summary(out2)
# plot(out2$MC) ##ploting Markov chain
# plot(out2)    ##ploting volatilities
 
### DCC-GARCH(1,1) with SSGED innovations ###
# out = bayesDccGarch(DaxCacNik, nSim=20000, tail_ini=1.5, errorDist=3)
# out2 = increaseSim(out, 10000) # more 10000 simulations
# out3 = window(out2, start=5000)  # removing the first 5000 as burning
# summary(out3)
#
### Plotting volatilities #
# plot(out3)
#
### Plotting Markov Chain #
# plot(out3$MC)
#
############################################################

######  log-likelihood function of GARCH(1,1) model with SST innovations ####
logLikDccGarch(Dax, omega=0.03, alpha=0.03, beta=0.8, gamma=0.7)$value

######  log-likelihood function of DCC-GARCH(1,1) model with SST innovations ####
logLikDccGarch(DaxCacNik, omega=c(0.02,0.01,0.03), alpha=c(0.02,0.01,0.03))$value

}
