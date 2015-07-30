\name{bayesDccGarch}
\alias{bayesDccGarch}
\alias{increaseSim}

\title{
 Bayesian Estimation of the DCC-GARCH(1,1) Model.
}
\description{
 Performs a Markov Chain for all parameters of the DCC-GARCH(1,1) Model.
}

\usage{
bayesDccGarch(mY, nSim = 10000, tail_ini = 8, omega_ini = rep(0.03, ncol(mY)), 
	alpha_ini = rep(0.03, ncol(mY)), beta_ini = rep(0.8, ncol(mY)), 
	a_ini = 0.03, b_ini = 0.8, gamma_ini = rep(1, ncol(mY)), 
	errorDist = 2, control = list())
	
increaseSim(outList, nSim=10000)	
}

\arguments{
  \item{mY}{a matrix of the data (\eqn{n \times k}).}
  \item{nSim}{length of Markov chain. Default: \code{10000}.}
  \item{tail_ini}{initial value of \eqn{\nu} parameter if \code{errorDist = 2} or initial value of \eqn{\delta} parameter if \code{errorDist = 3}. If \code{errorDist = 1} this arguments is not used.}
  \item{omega_ini}{a numeric vector (\eqn{k \times 1}) with the initial values of \eqn{\omega_i} parameters. Default: \code{rep(0.03, ncol(mY))}.}
  \item{alpha_ini}{a numeric vector (\eqn{k \times 1}) with the initial values of \eqn{\alpha_i} parameters. Default: \code{rep(0.03, ncol(mY))}.}
  \item{beta_ini}{a numeric vector (\eqn{k \times 1}) with the initial values of \eqn{\beta_i} parameters. Default: \code{rep(0.8, ncol(mY))}.}
  \item{a_ini}{a numeric value of the initial values of \eqn{a} parameter. Default: \code{0.03}.}
  \item{b_ini}{a numeric value of the initial values of \eqn{b} parameter. Default: \code{0.8}.}
  \item{gamma_ini}{a numeric vector (\eqn{k \times 1}) with the initial values of \eqn{\gamma_i} parameters. Default: \code{rep(1.0, ncol(mY))}.}
  \item{errorDist}{a probability distribution for errors. Use \code{errorDist=1} for \eqn{SSNorm}, \code{errorDist=2} for \eqn{SST} or \code{errorDist=3} for \eqn{SSGED}. Default: \code{2}.}
  \item{control}{list of control arguments (See *Details*).}
  \item{outList}{an output element of \code{bayesDccGarch} function.}
}

\details{
The \code{bayesDccGarch} performs a Markov Chain for all parameters of the model DCC-GARCH(1,1) (or GARCH(1,1) in the univariate case). 
There are three options of probability distributions for the error
component. These are the standardized skew versions of normal, t-student and ged distributions. See Fioruci et al (2014a) and Fioruci et al (2014b) for any detail.
 
The \code{control} argument can be used for define the prior hyper-parameters and the simulation algorithm parameters.  
It is a list that can supply any of the following components:
  \describe{	
  \item{$mu_tail}{the value of hyper-parameter \eqn{\mu_\nu} if \code{errorDist=2} or the hyper-parameter \eqn{\mu_\delta} if \code{errorDist=3}. Default: 8}
  \item{$mu_gamma}{a vector with the hyper-parameters \eqn{\mu_{\gamma_i}}. Default: \code{rep(0,ncol(mY)}}
  \item{$mu_omega}{a vector with the hyper-parameters \eqn{\mu_{\omega_i}}. Default: \code{rep(0,ncol(mY)}}
  \item{$mu_alpha}{a vector with the hyper-parameters \eqn{\mu_{\alpha_i}}. Default: \code{rep(0,ncol(mY)}}
  \item{$mu_beta}{a vector with the hyper-parameters \eqn{\mu_{\beta_i}}. Default: \code{rep(0,ncol(mY)}}
  \item{$mu_a}{the value of the hyper-parameter \eqn{\mu_a}. Default: 0}
  \item{$mu_b}{the value of the hyper-parameter \eqn{\mu_b}. Default: 0}
  \item{$sigma_tail}{the value of hyper-parameter \eqn{\sigma_\nu} if \code{errorDist=2} or the hyper-parameter \eqn{\sigma_\delta} if \code{errorDist=3}. Default: 10}
  \item{$sigma_gamma}{a vector with the hyper-parameters \eqn{\sigma_{\gamma_i}}. Default: \code{rep(1.25,ncol(mY)}}
  \item{$sigma_omega}{a vector with the hyper-parameters \eqn{\sigma_{\omega_i}}. Default: \code{rep(10,ncol(mY)}}
  \item{$sigma_alpha}{a vector with the hyper-parameters \eqn{\sigma_{\alpha_i}}. Default: \code{rep(10,ncol(mY)}}
  \item{$sigma_beta}{a vector with the hyper-parameters \eqn{\sigma_{\beta_i}}. Default: \code{rep(10,ncol(mY)}}
  \item{$sigma_a}{the value of the hyper-parameter \eqn{\sigma_a}. Default: 10}
  \item{$sigma_b}{the value of the hyper-parameter \eqn{\sigma_b}. Default: 10}
  
  \item{$simAlg}{the random walk Metropolis-Hasting algorithm update. Use \code{1} for update all parameters as one block, 
  use \code{2} for update one parameter for each time and use \code{3} for an automatic choice.}
  \item{$cholCov}{the cholesky decomposition matrix of the covariance matrix for simulation by one-block Metropolis-Hasting. It must to be passed if \code{control$simAlg=1}.}
  \item{$sdSim}{a vector with the standard deviations for simulation by one-dimensional Metropolis-Hasting. It must to be passed if \code{control$simAlg=2}.}
  
  \item{$print}{a logical variable for if the function should report the number of interactions in each 100 interactions or not. Default: TRUE} 
  }
}

\value{
A list with elements:
\item{$control }{ a list with the used \code{control} argument. }
\item{$MC }{ a element of \code{mcmc} class with the Markov Chain simulation of all parameters. (\R package \pkg{coda})}
\item{$H }{ a matrix with the Bayesian estimates of volatilities and co-volatilities.}
\item{$IC }{ the Bayesian estimate of Akaike Information Criterion, Bayesian Information Criterion and Deviance Information Criterion.}
\item{$elapsedTime }{ an object of class \code{proc_time} which is a numeric vector of length 5, containing the user, system, and total elapsed times of the process.}
}

\references{
Fioruci, J.A., Ehlers, R.S., Andrade Filho, M.G. \emph{Bayesian multivariate GARCH models with dynamic correlations and asymmetric error distributions}, Journal of Applied Statistics, 41(2), 320--331, 2014a. \url{http://dx.doi.org/10.1080/02664763.2013.839635}.

Fioruci, J.A., Ehlers, R.S., Louzada, F. \emph{BayesDccGarch - An Implementation of Multivariate GARCH DCC Models}, ArXiv e-prints, 2014b. \url{http://adsabs.harvard.edu/abs/2014arXiv1412.2967F}.
}

\author{Jose Augusto Fioruci, Ricardo Sandes Ehlers and Francisco Louzada}

\seealso{ \code{\link{bayesDccGarch-package}}, \code{\link{logLikDccGarch}}, \code{\link{plotVol}} }

\examples{

data(DaxCacNik)

mY = DaxCacNik[1:10,] # more data is necessary
Dax = DaxCacNik[,1]

### DCC-GARCH(1,1) ###
out = bayesDccGarch(mY, nSim=1000)
summary(out$MC)
plotVol(mY, out$H[,c("H_1,1","H_2,2","H_3,3")], c("DAX","CAC40","NIKKEI")) 


### GARCH(1,1) ###
# out = bayesDccGarch(Dax)
# summary(out$MC)
# plotVol(Dax, out$H, "DAX") 	


### DCC-GARCH(1,1) with SSGED innovations ###
## out = bayesDccGarch(mY, nSim=5000, tail_ini=1.5, errorDist=3)


### Informative Prior for alpha and beta parameters ###
# out = bayesDccGarch(Dax, alpha_ini=0.01, beta_ini=0.90, 
#	control=list(mu_alpha=0.01, sigma_alpha=0.001, mu_beta=0.90, sigma_beta=0.001))
# summary(out$MC)
# plot(out$MC)
# plotVol(Dax, out$H, "DAX") 	

### more simulations ###
# out1 = bayesDccGarch(mY)
# out2 = increaseSim(out1, 10000)
# summary(out2$MC)
# plot(out2$MC)
# plotVol(mY, out2$H[,c("H_1,1","H_2,2","H_3,3")], c("DAX","CAC40","NIKKEI"))

}

\keyword{ts}
\keyword{multivariate}
\keyword{dcc-garch}
\keyword{garch}
\keyword{models}

