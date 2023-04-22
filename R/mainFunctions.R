
memoryAllocation <- function(mY, errorDist, mu_omega, sigma_omega, mu_alpha, sigma_alpha, mu_beta, sigma_beta, mu_a, sigma_a, mu_b, sigma_b, mu_gamma, sigma_gamma, mu_tail, sigma_tail, print)
{
  out = .C("memoryAllocation", vY=as.numeric(t(mY)), vn=as.integer(nrow(mY)), vk=as.integer(ncol(mY)), verrorDist = as.integer(errorDist), vmu_omega = as.double(mu_omega), vsigma_omega = as.double(sigma_omega), vmu_alpha = as.double(mu_alpha), vsigma_alpha = as.double(sigma_alpha), vmu_beta = as.double(mu_beta), vsigma_beta = as.double(sigma_beta), vmu_a = as.double(mu_a), vsigma_a = as.double(sigma_a), vmu_b = as.double(mu_b), vsigma_b = as.double(sigma_b), vmu_gamma = as.double(mu_gamma), vsigma_gamma = as.double(sigma_gamma), vmu_tail = as.double(mu_tail), vsigma_tail = as.double(sigma_tail), vprint=as.integer(print) )
}

memoryDeallocation <- function(){
  out = .C("memoryDeallocation")
}

logPosterior_phi <- function(phi){
  if(length(phi)==5) phi = c(phi, -2, 1)

  .C("logPosterior_phi_R", phi = as.double(phi), value = numeric(1) )$value[1]
}

logLikDccGarch <- function(mY, omega=rep(0.03, ncol(mY)), alpha=rep(0.03, ncol(mY)), beta=rep(0.80, ncol(mY)), a=0.03, b=0.80, gamma=rep(1.00, ncol(mY)), tail=10, errorDist=2){
  if (missing(mY)) stop ("'mY' is missing")
  if (any(is.na(mY))) stop ("'mY' contains 'NA' values")
  mY=as.matrix(mY)
  if(length(omega)!=ncol(mY)) stop("The omega argument must be a numeric vector and satisfy the condition \"length(omega)=ncol(mY)\".")
  if(length(alpha)!=ncol(mY)) stop("The alpha argument must be a numeric vector and satisfy the condition \"length(alpha)=ncol(mY)\".")
  if(length(beta)!=ncol(mY)) stop("The beta argument must be a numeric vector and satisfy the condition \"length(beta)=ncol(mY)\".")
  if(length(gamma)!=ncol(mY)) stop("The gamma argument must be a numeric vector and satisfy the condition \"length(gamma)=ncol(mY)\".")
  if(length(a)!=1) stop("The \"a\" argument must satisfy the condition \"length(a)=1\".")
  if(length(b)!=1) stop("The \"b\" argument must satisfy the condition \"length(b)=1\".")
  if(!any(errorDist==1:3)) stop("The errorDist argument must be one element of the numbers 1, 2 and 3.")
  if(any(c(omega,alpha,beta,a,b,gamma)<0) || any(c(alpha+beta,a+b)>1)) stop("There is at least one parameter out of parametric interval.")

  out = .C("logLikelihood_R", y=as.numeric(t(mY)), omega=as.double(omega), alpha=as.double(alpha), beta=as.double(beta), a=as.double(a), b=as.double(b), gamma=as.double(gamma), tail=as.double(tail), verrorDist=as.integer(errorDist), vn=as.integer(nrow(mY)), vk=as.integer(ncol(mY)), value=numeric(1), vmH=numeric(nrow(mY)*ncol(mY)*ncol(mY)))

  out2= list()
  out2$H = matrix(out$vmH, ncol=ncol(mY)^2, byrow=T)
  out2$value = out$value

  return(out2)
}


updateBayesDccGarch <- function(x, n_mcmc_update=Inf){
  mY = x$control$data
  n = nrow(mY)
  k = ncol(mY)
  n_mcmc = nrow(x$MC)

  if(k==1){
	x$MC = cbind(x$MC, 0.1, 0.1) ## include a,b parameters
  }

  if(n_mcmc > n_mcmc_update){
    ss =  round(seq.int(from=1, to=n_mcmc, length.out=n_mcmc_update)) #sample.int(n=n_mcmc, size = 2000, replace = FALSE)
  }else{
	  ss = 1:n_mcmc
  }

  out = .C("loglike_matrix_R",
    vY = as.numeric(t(mY)),
    vn = as.integer(nrow(mY)),
    vk = as.integer(ncol(mY)),
    vmcmc = as.numeric(t(x$MC[ss,])),
    n_mcmc = as.integer(length(ss)), #as.integer(n_mcmc),
    verrorDist = as.integer(x$control$errorDist),
    meanLogLike = numeric(1),
    vmMeanH = numeric(nrow(mY)*ncol(mY)*ncol(mY))
  )

  if(k==1){
	  x$MC = as.mcmc( x$MC[,1:(4*k+1)] ) ## exclude a,b parameters
  }

  D_mean = -2 * out$meanLogLike
  par_est = apply(X=x$MC, MARGIN=2, FUN=mean)
  tail_est = par_est[1]; gamma_est = par_est[2+ 4*(0:(k-1))]; omega_est = par_est[3+ 4*(0:(k-1))];
  alpha_est = par_est[4+ 4*(0:(k-1))]; beta_est = par_est[5+ 4*(0:(k-1))];
  if(k>1){ a_est = par_est[2+4*k]; b_est = par_est[3+4*k];  }else{ a_est = 0.1; b_est = 0.1;}

  D_at_mean = -2*logLikDccGarch(mY, omega_est, alpha_est, beta_est, a_est, b_est, gamma_est, tail_est, x$control$errorDist)$value

  if( k == 1 ){ npar = 5 }else{ npar = 4*k +3 }

  if(x$control$errorDist==1){
    EAIC = D_mean + 2*(npar-1)
	  EBIC = D_mean + (npar-1)*log(n)
  }else{
	  EAIC = D_mean + 2*npar
	  EBIC = D_mean + npar*log(n)
  }
  DIC = D_at_mean + 2*(D_mean - D_at_mean)

  IC = matrix(c(EAIC, EBIC, DIC), nrow=3, ncol=1)
  rownames(IC) = c("EAIC", "EBIC", "DIC")
  colnames(IC) = c("Estimative")
  x$IC = IC

  Hname = rep("NA", k*k)
  for(i in 1:k){ for(j in 1:k){ Hname[(i-1)*k + j] = paste("H_",i,",",j,sep="") } }

  #Hname = colnames(x$H)
  x$H = matrix(out$vmMeanH, ncol=ncol(mY)^2, byrow=T)
  colnames(x$H) = Hname

  ### Predicted Co-Volatilities and Correlation matrix for y_{n+1}
  H_n1 <- matrix( tail(logLikDccGarch(rbind(mY,rep(0,ncol(mY))), omega_est, alpha_est, beta_est, a_est, b_est, gamma_est, tail_est, x$control$errorDist)$H,1), k, k )
  colnames(H_n1) = rownames(H_n1) = colnames(mY)
  R_n1 <- H_n1

  ### Compute R matrix and R_n1
  DCC <- matrix(NA, nrow = nrow(x$H), ncol=ncol(x$H)) ### Dynamic Conditional Correlation
  colnames(DCC) = colnames(x$H)

  for(i in 1:k){
	col_i_name <- paste0("H_",i,",",i)

	 for(j in 1:k){
		col_j_name <- paste0("H_",j,",",j)

		col_ij_name <- paste0("H_",i,",",j)

		DCC[,col_ij_name] <- x$H[,col_ij_name] / (sqrt(x$H[,col_i_name]) * sqrt(x$H[,col_j_name]))

		R_n1[i,j] <- H_n1[i,j]/(sqrt(H_n1[i,i])*sqrt(H_n1[j,j]))
	 }
  }

  colnames(DCC) <- colnames(x$R)
  x$R <- DCC
  x$H_n1 <- H_n1
  x$R_n1 <- R_n1

  #x$logLike = out$meanLogLike

  return(x)
}


MH_oneDimension <- function(phi_ini, k, n_sim, sd_phi_sim){

  if( k == 1 )
    npar = 5
  else
    npar = 4*k +3

  vMC = .C("MH_oneDimension", phi = as.double(phi_ini), sd_phi_sim = as.double(sd_phi_sim), n_sim = as.integer(n_sim), vMC=numeric(n_sim*npar) )$vMC

  MC = as.mcmc( matrix(vMC, ncol=npar, byrow=T) )
  colnames(MC) = paste("phi_", 1:npar, sep="")

  return(MC)
}

MH_oneBlock <- function(phi_ini, k, n_sim, chol_cov_phi_sim){

  if( k == 1 )
    npar = 5
  else
    npar = 4*k +3

  vMC = .C("MH_oneBlock", phi = as.double(phi_ini), vChol_Cov_phi_sim = as.numeric(t(chol_cov_phi_sim)), n_sim = as.integer(n_sim), vMC=numeric(n_sim*npar) )$vMC

  MC = as.mcmc( matrix(vMC, ncol=npar, byrow=T) )
  colnames(MC) = paste("phi_", 1:npar, sep="")

  return(MC)
}


bayesDccGarch <- function(mY, nSim=10000, tail_ini=8, omega_ini=0.1*diag(var(mY)), alpha_ini=rep(0.05, ncol(mY)), beta_ini=rep(0.85, ncol(mY)), a_ini=0.04, b_ini=0.8, gamma_ini=rep(1.0, ncol(mY)), errorDist=2, control=list() )
{
	ptm <- proc.time()

	if (missing(mY)) stop ("'mY' is missing")
  if (any(is.na(mY))) stop ("'mY' contains 'NA' values")
	if (errorDist==3 && tail_ini==8){ tail_ini=2 }
  mY=as.matrix(mY)
  n = nrow(mY)
	k = ncol(mY)
	if( k == 1 ){ npar = 5 }else{ npar = 4*k +3 }

	Lout = try( logLikDccGarch(mY, omega_ini, alpha_ini, beta_ini, a_ini, b_ini, gamma_ini, tail_ini, errorDist) )

	if("try-error" %in% class(Lout)) stop("The likelihood function can not be computed for the initial values.")
	if(!is.list(control)) stop("control must be an object of list class.")
	if(is.null(control$data)){ control$data = mY }
	if(is.null(control$errorDist)){ control$errorDist = errorDist }
	if(is.null(control$mu_a)){ control$mu_a = 0	}else{ if( length(control$mu_a)!=1 ){ stop("control$mu_a must be a numeric vector with length equal to 1.") } }
  if(is.null(control$mu_b)){ control$mu_b = 0 }else{ if( length(control$mu_b)!=1 ){ stop("control$mu_b must be a numeric vector with length equal to 1.") } }
	if(is.null(control$mu_omega)){ control$mu_omega = rep(0,k) }else{ if( length(control$mu_omega)!=k ){ stop("control$mu_omega must be a numeric vector with length equal to 'ncol(mY)'.")} }
	if(is.null(control$mu_alpha)){ control$mu_alpha = rep(0,k) }else{ if( length(control$mu_alpha)!=k ){ stop("control$mu_alpha must be a numeric vector with length equal to 'ncol(mY)'.")} }
	if(is.null(control$mu_beta)){ control$mu_beta = rep(0,k) }else{ if( length(control$mu_beta)!=k ){ stop("control$mu_beta must be a numeric vector with length equal to 'ncol(mY)'.")} }
	if(is.null(control$mu_gamma)){ control$mu_gamma = rep(0,k) }else{ if( length(control$mu_gamma)!=k ){ stop("control$mu_gamma must be a numeric vector with length equal to 'ncol(mY)'.")} }
  if(is.null(control$mu_tail)){ control$mu_tail = 0 }else{ if( length(control$mu_tail)!=1 ){ stop("control$mu_tail must be a numeric vector with length equal to 1.")} }
	if(is.null(control$sigma_a)){ control$sigma_a = 10 }else{ if( length(control$sigma_a)!=1 ){ stop("control$sigma_a must be a numeric vector with length equal to 1.")} }
  if(is.null(control$sigma_b)){ control$sigma_b = 10 }else{ if( length(control$sigma_b)!=1 ){ stop("control$sigma_b must be a numeric vector with length equal to 1.")} }
	if(is.null(control$sigma_omega)){ control$sigma_omega = rep(10,k) }else{ if( length(control$sigma_omega)!=k ){ stop("control$sigma_omega must be a numeric vector with length equal to 'ncol(mY)'.")} }
	if(is.null(control$sigma_alpha)){ control$sigma_alpha = rep(10,k) }else{ if( length(control$sigma_alpha)!=k ){ stop("control$sigma_alpha must be a numeric vector with length equal to 'ncol(mY)'.")} }
	if(is.null(control$sigma_beta)){ control$sigma_beta = rep(10,k) }else{ if( length(control$sigma_beta)!=k ){ stop("control$sigma_beta must be a numeric vector with length equal to 'ncol(mY)'.")} }
	if(is.null(control$sigma_gamma)){ control$sigma_gamma = rep(10,k) }else{ if( length(control$sigma_gamma)!=k ){ stop("control$sigma_gamma must be a numeric vector with length equal to 'ncol(mY)'.")} }
  if(is.null(control$sigma_tail)){ control$sigma_tail = 1.25 }else{ if( length(control$sigma_tail)!=1 ){ stop("control$sigma_tail must be a numeric vector with length equal to 1.")} }
	if(is.null(control$print) ){ control$print = 1 }else{ if(!any(control$print == 0:1)) control$print = 1 }
	if(is.null(control$simAlg)){ control$simAlg = 3 }else{ if(!any(control$simAlg == 1:3)) control$simAlg = 3 }
	if(control$simAlg == 1){
		option = 1
		if( is.null(control$cholCov) ){ stop("control$cholCov was not found.")
		}else{ if(!is.matrix(control$cholCov)){ stop("control$cholCov must be a positive definite matrix of dimension 5x5 if k=1 or (4k+3)x(4k+3) if k>1, where 'k=ncol(mY)'")
		       }else{  if(ncol(control$cholCov)!=npar || nrow(control$cholCov)!=npar){ stop("control$cholCov must be a positive definite matrix of dimension 5x5 if k=1 or (4k+3)x(4k+3) if k>1, where 'k=ncol(mY)'")
						}
				}
		}
	}
	if(control$simAlg == 2 && !is.null(control$sdSim)){
		option = 2
		if(length(control$sdSim)!=npar) stop("control$sdSim must be a numeric positive vector with length 5 if k=1 or 4k+3 if k>1, where 'k=ncol(mY)'")
	}
	if(control$simAlg == 3){ if(k==1){control$sdSim = c(0.2, 0.2, 0.4, 0.4, 0.2) }else{control$sdSim = c(0.2, rep(c(0.2, 0.4, 0.4, 0.2), k), 1.0, 1.5)} }
	if(is.null(control$lambda)){control$lambda = 0.4}
	if(is.null(control$nPilotSim)){control$nPilotSim = 1000}



  phi_ini = .C("original_to_real_scale",
                  phi = numeric(4*k+3),
                  omega = as.double(omega_ini),
                  alpha = as.double(alpha_ini),
                  beta = as.double(beta_ini),
                  a = as.double(a_ini),
                  b = as.double(b_ini),
                  gamma = as.double(gamma_ini),
                  tail = as.double(tail_ini),
                  k = as.integer(k),
                  errorDist = as.integer(errorDist)
              )$phi


  memoryAllocation( mY, errorDist, control$mu_omega, control$sigma_omega, control$mu_alpha,
					control$sigma_alpha, control$mu_beta, control$sigma_beta,
					control$mu_a, control$sigma_a, control$mu_b, control$sigma_b,
					control$mu_gamma, control$sigma_gamma, control$mu_tail, control$sigma_tail,
					control$print )


    # logPosterior_phi(phi_ini)
    # .C("printGlobalMatrix")

	if( control$simAlg == 3 ){
		writeLines("Maximizing the log-posterior density function.")
		if(k==1) phi_ini = phi_ini[1:5]
		opt = optim(par = phi_ini, fn = logPosterior_phi, hessian=TRUE, control = list(fnscale=-abs(Lout$value))) #optim(par = phi_ini, fn = logPosterior_phi, hessian=TRUE, control = list(fnscale=-1))
		if(k==1) phi_ini = c(phi_ini,-2,1)

		if(opt$convergence==0 || opt$convergence==1){
			writeLines("Done.")
			phi_ini = opt$par
			if(k==1) phi_ini = c(phi_ini,-2,1)

			control$cholCov = try( t( chol( solve( -opt$hessian ) ) ) , silent = TRUE )
			# if( class(control$cholCov) == "try-error"){
				# nearPD_Hessian = nearPD( -opt$hessian )

				# if( nearPD_Hessian$converged ){
				  # control$cholCov = try( t( chol( solve( nearPD_Hessian$mat ) ) ) , silent = TRUE )
				# }
			# }

			if( !("try-error" %in% class(control$cholCov)) ){ #if( class(control$cholCov) != "try-error"){
				option = 1
			}else{
				writeLines("One approximation for covariance matrix of parameters cannot be directly computed through the hessian matrix.")
				option = 2
			}
		}else{
			writeLines("The optim function convergence was not obtained.")
			option = 2
		}
	}

    repeat{

        if(option==1 && control$simAlg!=2){
            #writeLines("\nStarting the first simulation option")

			#if(control$simAlg==3){
				writeLines("Calibrating the Lambda coefficient:")

				cont=1
				proceed=TRUE
				while(proceed){
					writeLines(paste("lambda:", control$lambda))
					MC_phi = MH_oneBlock(phi_ini, k, n_sim=control$nPilotSim, control$lambda*control$cholCov)
					acceptRate = 1-rejectionRate( MC_phi )[5]
					writeLines(paste("Accept Rate:", round(acceptRate,2)))
					if(acceptRate>.20 && acceptRate<.50 || cont>5){
						proceed=FALSE
					}else{
						cont=cont+1
						if(acceptRate < .20){ control$lambda = control$lambda*0.8  }
						if(acceptRate > .50){ control$lambda = control$lambda*1.2 }
					}
				}

				writeLines("Done.")
            #}

            writeLines("Starting the simulation by one-block random walk Metropolis-Hasting algorithm.")
            MC_phi = MH_oneBlock(phi_ini, k, n_sim=nSim, control$lambda*control$cholCov)
            writeLines("Done.")
            break
        }

      if(option==2 && control$simAlg!=1 ){
            #writeLines("\nStarting the second simulation option.")

			if(control$simAlg == 3){
				writeLines("Calibrating the standard deviations for simulation:")

				proceed=TRUE
				cont=1
				while(proceed){
					MC_phi = MH_oneDimension(phi_ini, k=k, n_sim=control$nPilotSim, control$sdSim)
					acceptRate = 1-rejectionRate( MC_phi )
					writeLines("Accept Rate:")
					print(round(acceptRate,2))
					if(errorDist == 1){
						acceptRate[1] = .30
					}
					lowerRate = which(acceptRate<.15)
					upperRate = which(acceptRate>.50)
					if(length(lowerRate) == 0 && length(upperRate) == 0 || cont>5){
						proceed=FALSE
					}else{
						cont=cont+1
						control$sdSim[lowerRate] = control$sdSim[lowerRate]*0.6
						control$sdSim[upperRate] = control$sdSim[upperRate]*1.4
					}
				}

				if(errorDist == 1){
					MC_phi[,1] = rnorm(control$nPilotSim)
				}

				writeLines("Computing the covariance matrix of pilot sample.")
				control$cholCov = try( t( chol( cov(MC_phi) ) ) , silent = TRUE )
				if( !("try-error" %in% class(control$cholCov)) ){  #if( class(control$cholCov) != "try-error"){
					writeLines("Done.")
					option = 1
				}else{
					control$cholCov = NULL
					writeLines("The approximately covariance matrix is not positive definitely.")
					option = 2
				}

      }

			if(option == 2){
				writeLines("Starting the simulation by one-dimensional random walk Metropolis-Hasting algorithm.")
				MC_phi = MH_oneDimension(phi_ini, k=k, n_sim=nSim, sd_phi_sim=control$sdSim)
				writeLines("Done.")
				break
			}

        }

    }

	vMMeanH = .C("getMeanH", vMMeanH= numeric(n*k*k) )$vMMeanH
	mMeanH = matrix(vMMeanH,nrow=n,ncol=k*k,byrow=TRUE)

	name = name2 =rep("NA", k*k)
	for(i in 1:k){
		for(j in 1:k){
			name[(i-1)*k + j] = paste("H_",i,",",j,sep="")
			name2[(i-1)*k + j] = paste("R_",i,",",j,sep="")
		}
	}
	colnames(mMeanH) = name

	logLike_mean = .C("getLogLikelihood_mean", value=numeric(1))$value


	memoryDeallocation()

	### Parameters inverse transformation
	name = rep("NA", npar)
    MC = matrix(NA, nrow=nSim, ncol=npar)
    if( errorDist==1 || errorDist==3){
        MC[,1] = exp(MC_phi[,1])
        name[1] = "delta"
    }
    if( errorDist==2){
        MC[,1] = exp(MC_phi[,1])+2
        name[1] = "nu"
    }
    for(i in 1:k){
        MC[,4*(i-1)+2:5] = cbind( exp(MC_phi[,4*(i-1) +2:3]), exp(MC_phi[,4*(i-1) +4:5])/(1+exp(MC_phi[,4*(i-1) +4:5])) )
        name[4*(i-1)+2:5] = paste( c("gamma", "omega", "alpha", "beta"), "_", i, sep="")
    }
    if(k>1){
        MC[,4*k +2:3] =  exp(MC_phi[,4*k +2:3])/(1+exp(MC_phi[,4*k +2:3]))
        name[4*k +2:3] = c("a", "b")
    }
    colnames(MC) = name

	### Compute criteria informations
	D_mean = -2*logLike_mean
	par_est = apply(X=MC, MARGIN=2, FUN=mean)
	tail_est = par_est[1]; gamma_est = par_est[2+ 4*(0:(k-1))]; omega_est = par_est[3+ 4*(0:(k-1))];
	alpha_est = par_est[4+ 4*(0:(k-1))]; beta_est = par_est[5+ 4*(0:(k-1))];
	if(k>1){ a_est = par_est[2+4*k]; b_est = par_est[3+4*k];  }else{ a_est = 0.1; b_est = 0.1;}

	D_at_mean = -2*logLikDccGarch(mY, omega_est, alpha_est, beta_est, a_est, b_est, gamma_est, tail_est, errorDist)$value

	if(errorDist==1){
		EAIC = D_mean + 2*(npar-1)
		EBIC = D_mean + (npar-1)*log(n)
	}else{
		EAIC = D_mean + 2*npar
		EBIC = D_mean + npar*log(n)
	}
	DIC = D_at_mean + 2*(D_mean - D_at_mean)

	IC = matrix(c(EAIC, EBIC, DIC), nrow=3, ncol=1)
	rownames(IC) = c("EAIC", "EBIC", "DIC")
	colnames(IC) = c("Estimative")

	### Predicted Co-Volatilities and Correlation matrix for y_{n+1}
	H_n1 <- matrix( tail(logLikDccGarch(rbind(mY,rep(0,ncol(mY))), omega_est, alpha_est, beta_est, a_est, b_est, gamma_est, tail_est, errorDist)$H,1), k, k )
	colnames(H_n1) = rownames(H_n1) = colnames(mY)
	R_n1 <- H_n1

	### Compute R matrix and R_n1
	DCC <- matrix(NA, nrow = nrow(mMeanH), ncol=ncol(mMeanH)) ### Dynamic Conditional Correlation
	colnames(DCC) = colnames(mMeanH)

	for(i in 1:k){
		col_i_name <- paste0("H_",i,",",i)

		for(j in 1:k){
			col_j_name <- paste0("H_",j,",",j)

			col_ij_name <- paste0("H_",i,",",j)

			DCC[,col_ij_name] <- mMeanH[,col_ij_name] / (sqrt(mMeanH[,col_i_name]) * sqrt(mMeanH[,col_j_name]))

			R_n1[i,j] <- H_n1[i,j]/(sqrt(H_n1[i,i])*sqrt(H_n1[j,j]))
		}
	}
	colnames(DCC) = name2

  out = list()
	out$control = control
  out$MC = as.mcmc(MC)
	out$MC_phi = as.mcmc(MC_phi)
	out$H = mMeanH
	out$R = DCC
	out$H_n1 = H_n1
	out$R_n1 = R_n1
	out$IC = IC
	# out$logLike = logLike_mean
  out$elapsedTime = proc.time() - ptm

  return(structure(out,class="bayesDccGarch"))
}


increaseSim <- function(x, nSim=10000)
{
  if( !("bayesDccGarch" %in% class(x)) ){ stop("Error: argument x is not one element of 'bayesDccGarch' class") }
  control = x$control
  k = ncol(control$data)
  nMC = nrow(x$MC)
  tail_ini = x$MC[nMC,1]
  gamma_ini = x$MC[nMC,2 + 4*(0:(k-1))]
  omega_ini = x$MC[nMC,3 + 4*(0:(k-1))]
  alpha_ini = x$MC[nMC,4 + 4*(0:(k-1))]
  beta_ini =  x$MC[nMC,5 + 4*(0:(k-1))]
  if(k > 1){
    a_ini = x$MC[nMC, 2 + 4*k]
    b_ini = x$MC[nMC, 3 + 4*k]
  }else{
    a_ini = 0.03
	b_ini = 0.8
  }

  if(!is.null(control$cholCov)){ control$simAlg = 1}else{ control$simAlg = 2 }

  outNewSim = bayesDccGarch(control$data, nSim, tail_ini, omega_ini, alpha_ini, beta_ini, a_ini, b_ini, gamma_ini, control$errorDist, control)

  x$H = (nMC*x$H + nSim*outNewSim$H)/(nMC + nSim)
  x$R = (nMC*x$R + nSim*outNewSim$R)/(nMC + nSim)
  x$H_n1 = (nMC*x$H_n1 + nSim*outNewSim$H_n1)/(nMC + nSim)
  x$R_n1 = (nMC*x$R_n1 + nSim*outNewSim$R_n1)/(nMC + nSim)
  x$IC = (nMC*x$IC + nSim*outNewSim$IC)/(nMC + nSim)
  x$MC = as.mcmc( rbind(x$MC, outNewSim$MC) )
  x$MC_phi = as.mcmc( rbind(x$MC_phi, outNewSim$MC_phi) )
  x$elapsedTime = x$elapsedTime + outNewSim$elapsedTime

  return(x)
}

window.bayesDccGarch <- function(x, start = NULL, end = NULL, thin = NULL,...){
	if( !("bayesDccGarch" %in% class(x)) ){ stop("Error: argument x is not a element of 'bayesDccGarch' class") }
	nMC = nrow(x$MC)
	if(is.null(start)){start = 1}
	if(is.null(end)){end = nMC}
	if(is.null(thin)){thin = 1}

	x$MC = window(as.mcmc(x$MC), start = start, end = end, thin = thin,...)
	x$MC_phi = window(as.mcmc(x$MC_phi), start = start, end = end, thin = thin,...)

	return( updateBayesDccGarch(x) )
}


plotVol <- function(mY, vol, ts.names=paste("TS_", 1:ncol(mY), sep=""), colors = c("grey","red"), ...){
  if (missing(mY)) stop ("'mY' is missing.")
  if (missing(vol)) stop ("'vol' is missing.")
  if (any(is.na(mY))) stop ("'mY' contains 'NA' values.")
  if (any(is.na(vol))) stop ("'vol' contains 'NA' values.")
  mY=as.matrix(mY)
  vol=as.matrix(vol)
  if(!all(dim(mY) == dim(vol))) stop ("the dimensions of 'mY' and 'vol' arguments are not the same.")
  if (ncol(mY) != length(ts.names)) stop ("the number of columns of 'ts.name' is different of length of 'ts.name'.")

  n = nrow(mY)
  k = ncol(mY)

  oldpar <- par(no.readonly = TRUE)    # code line i
  on.exit(par(oldpar))            # code line i + 1

  par(mfrow=c(k,1), mar= c(4, 4, 2, 2))
  for(i in 1:k){
	ylim = c(0, c(max(c(abs(mY[,i]), vol[,i]) )) )
    ts.plot(abs(mY[,i]),type='h',ylab=paste("|returns| of ",ts.names[i], sep=""), ylim=ylim, col=colors[1], ...)
    lines(vol[,i],col=colors[2]) #lines(exp(vol[,i]/2),col="red")
  }
}


summary.bayesDccGarch <- function(object,...){
	summary_mcmc = summary(as.mcmc(object$MC), ...)
	return(summary_mcmc)
}

print.bayesDccGarch <- function(x,...)
{
	dist_names = c("Standard Skew Normal", "Standard Skew t-Student", "Standard Skew GED")

	cat("Bayesian DCC-GARCH(1,1) model\n")
	cat( paste("Error Distribution:", dist_names[x$control$errorDist], "\n" ) )

	print( summary(x, ...) )

	cat( "3. Information criterions:\n" )

	print( x$IC )
	# summary.out = summary(x, ...)

	# cat("Statistics:\n")
	# print( summary.out$statistics )
	# cat("\nQuantiles:\n")
	# print( summary.out$quantiles )
}

plot.bayesDccGarch <- function(x, ts.names=NULL, colors = c("grey","red"), ...){

  mY = as.matrix( x$control$data )
  nTs = ncol(mY)
  vol = as.matrix(x$H[,paste("H_", 1:nTs, ",", 1:nTs, sep="")])
  if( is.null(ts.names) ){
  	if( is.null(colnames(mY)) ){
  	  ts.names = paste("TS_", 1:nTs, sep="")
  	}else{
  	  ts.names = colnames(mY)
  	}
  }

  plotVol(mY, vol, ts.names, colors, ...)
}




## included in version 2.3
update.bayesDccGarch <- function(object,..., mY_new){
  x = object
  if( !("bayesDccGarch" %in% class(x)) ){ stop("Error: argument x is not a element of 'bayesDccGarch' class") }
  mY_new = as.matrix(mY_new)
  mY = rbind(x$control$data, mY_new)
  x$control$data = mY
  return( updateBayesDccGarch(x) )
}



## included in version 3.0
#' Bayesian forecast for volatilities and coditional correlations
#' @aliases predict
#' @param object a bayesDccGarch object
#' @param ... default argument of predict function, not used
#' @param n_ahead number of steps ahead forecast
#' @param bayes a boolean. If True, then the forecast is calculated as being the average of the forecasts across all states in the Markov chain (much slower). If False then predictions are calculated using estimation parameters (much faster).
#' @return A list with elements \code{H} and \code{R}
#' @references
#' Engle, R.F. and Sheppard, K. Theoretical and empirical properties of dynamic conditional correlation multivariate GARCH, 2001, NBER Working Paper.
#' @examples
#' \donttest{
#' out = bayesDccGarch(DaxCacNik)
#' predict.bayesDccGarch(out, n_ahead=5)
#' }
predict.bayesDccGarch <- function(object,..., n_ahead=5, bayes=T){
  x = object
  if( !("bayesDccGarch" %in% class(x)) ){ stop("Error: argument x is not a element of 'bayesDccGarch' class") }
  mY = x$control$data
  n = nrow(mY)
  k = ncol(mY)

  R = as.numeric( cor( x$control$data ) ) # unconditional correlation matrix

  R_n1 = as.numeric(x$R_n1) # forecast one step ahead
  H_n1 = as.numeric(x$H_n1) # forecast one step ahead
  names(H_n1) = colnames(x$H)

  R_forec_out = H_forec_out = matrix(0, nrow=n_ahead, ncol=length(R_n1))
  colnames(R_forec_out) = colnames(x$R)
  colnames(H_forec_out) = colnames(x$H)
  rownames(R_forec_out) = rownames(H_forec_out) = paste0(1:n_ahead," step ahead")

  R_forec = H_forec = matrix(NA, nrow=n_ahead, ncol=length(R_n1))
  colnames(R_forec) = colnames(x$R)
  colnames(H_forec) = colnames(x$H)

  nSim = nrow(x$MC)
  ii=1
  repeat{

    idx = ii
    if(!bayes){
      idx = 1:nSim
    }

    if(k > 1){
      omega = colMeans( x$MC[idx, 3 + 4*(0:(k-1)), drop=F] )
      alpha = colMeans( x$MC[idx, 4 + 4*(0:(k-1)), drop=F] )
      beta =  colMeans( x$MC[idx, 5 + 4*(0:(k-1)), drop=F] )
      a = mean(x$MC[idx, 'a', drop=F])
      b = mean(x$MC[idx, 'b', drop=F])
    }else{
      omega = mean( x$MC[idx, 3, drop=F] )
      alpha = mean( x$MC[idx, 4, drop=F] )
      beta =  mean( x$MC[idx, 5, drop=F] )
      a = b = 0
    }

    R_forec[1,] = R_n1
    H_forec[1,] = H_n1
    if(n_ahead > 1){
      for(i in 1:length(R_n1)){
        R_forec[,i] = ( (a+b)^(0:(n_ahead-1)) )*( R_n1[i] - R[i] ) + R[i]
      }

      for(i in 1:k){
        H_forec[ ,paste0("H_",i,",",i)] = omega[i]*( (alpha[i]+beta[i])^(0:(n_ahead-1)) -1)/(alpha[i]+beta[i]-1)  + ((alpha[i]+beta[i])^(0:(n_ahead-1)))*H_n1[paste0("H_",i,",",i)]
      }
    }

    if(k > 1){
      for(i in 1:(k-1)){
        for(j in (i+1):k){
          H_forec[,paste0("H_",i,",",j)] = R_forec[,paste0("R_",i,",",j)]*sqrt(H_forec[,paste0("H_",i,",",i)])*sqrt(H_forec[,paste0("H_",j,",",j)])
          H_forec[,paste0("H_",j,",",i)] = H_forec[,paste0("H_",i,",",j)]
        }
      }
    }

    R_forec_out = R_forec_out + R_forec/ifelse(bayes, nSim, 1)
    H_forec_out = H_forec_out + H_forec/ifelse(bayes, nSim, 1)

    ii = ii+1
    if( !bayes || ii > nSim ){
      break
    }

  }

  return( list(H=H_forec_out, R=R_forec_out) )
}


