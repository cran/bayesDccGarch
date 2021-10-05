// RCMD SHLIB bayesDccGarch.c

#include<stdio.h>
#include<stdlib.h>
#include<R.h>
#include<Rmath.h>
#include<limits.h>
#include<math.h>
#include<time.h>
#include "matrixOperations.h"


double **y; //data matrix
int n, k;
int errorDist;
int print;
double **mH; //matrix of Ht's for t=1,\dots,n
double **mMeanH; // matrix to save the estimative of H_t's for t=1,...,n
//double **mcholH1; //matrix of Ht^{-1/2}'s for t=1,\dots,n
double **MEs; //Matrix Errors
mat R; // The unconditional covariance matrix
mat H, cholH, cholH1; // H_t, H_t^{-1} and H_t^{-1/2} matrix
mat Q; // Q_t matrix
vec mu_omega;
vec sigma_omega;
vec mu_alpha;
vec sigma_alpha;
vec mu_beta;
vec sigma_beta;
vec mu_a;
vec sigma_a;
vec mu_b;
vec sigma_b;
vec mu_gamma;
vec sigma_gamma;
vec mu_tail;
vec sigma_tail;
double logLikelihood_last, logLikelihood_old, logLikelihood_mean; //used for compute the DIC information

void getMeanH( double *vMMeanH ){
	matrix_to_vector(n, k*k, mMeanH, vMMeanH);
}

void getLogLikelihood_mean(double *value){
	value[0] = logLikelihood_mean;
}

void printGlobalMatrix(){
	Rprintf("\n\ny:\n");
	RprintMatrix(n, k, y);

	Rprintf("\n\nR:\n");
	RprintMatrix(k, k, R);

	Rprintf("\n\nH1:\n");
	RprintMatrix(k, k, H);

	Rprintf("\n\nQ1:\n");
	RprintMatrix(k, k, Q);

	Rprintf("\n\ncholH:\n");
	RprintMatrix(k, k, cholH);

	Rprintf("\n\ncholH1:\n");
	RprintMatrix(k, k, cholH1);

	Rprintf("\n\nmu_tail:\n");
	RprintVector(1, mu_tail);

	Rprintf("\n\nmu_gamma:\n");
	RprintVector(k, mu_gamma);

	Rprintf("\n\nmu_omega:\n");
	RprintVector(k, mu_omega);

	Rprintf("\n\nmu_alpha:\n");
	RprintVector(k, mu_alpha);

	Rprintf("\n\nmu_beta:\n");
	RprintVector(k, mu_beta);

	Rprintf("\n\nmu_a:\n");
	RprintVector(1, mu_a);

	Rprintf("\n\nmu_b:\n");
	RprintVector(1, mu_b);

	Rprintf("\n\nsigma_tail:\n");
	RprintVector(1, sigma_tail);

	Rprintf("\n\nsigma_gamma:\n");
	RprintVector(k, sigma_gamma);

	Rprintf("\n\nsigma_omega:\n");
	RprintVector(k, sigma_omega);

	Rprintf("\n\nsigma_alpha:\n");
	RprintVector(k, sigma_alpha);

	Rprintf("\n\nsigma_beta:\n");
	RprintVector(k, sigma_beta);

	Rprintf("\n\nsigma_a:\n");
	RprintVector(1, sigma_a);

	Rprintf("\n\nsigma_b:\n");
	RprintVector(1, sigma_b);
}

void zeroGlobalMatrix(){
	mat_zero(H, k);
	mat_zero(cholH, k);
	mat_zero(cholH1, k);
	mat_zero(Q, k);
	mat_zero(R, k);
}

void memoryAllocation(double *vY,
					int *vn,
					int *vk,
					int *verrorDist,
					double *vmu_omega,
					double *vsigma_omega,
					double *vmu_alpha,
					double *vsigma_alpha,
					double *vmu_beta,
					double *vsigma_beta,
					double *vmu_a,
					double *vsigma_a,
					double *vmu_b,
					double *vsigma_b,
					double *vmu_gamma,
					double *vsigma_gamma,
					double *vmu_tail,
					double *vsigma_tail,
					int *vprint
					){

	print = vprint[0];
	n = vn[0];
	k = vk[0];
	errorDist = verrorDist[0];

	H=mat_new(k); Q=mat_new(k); R=mat_new(k);
	cholH =mat_new(k);
	cholH1=mat_new(k); //Cholesky decomposition of H1
	matrix_new(k, n, &MEs); //allocation of errors matrix (Note that the order is k X n )

	matrix_new(n, k, &y); // allocation of the data matrix
	vector_to_matriz(n, k, vY, y);

	matrix_new(n, k*k, &mH); // matrix to save the H_t's for t=1,...,n
	matrix_new(n, k*k, &mMeanH); // matrix to save the estimative of H_t's for t=1,...,n

	mu_tail = vec_copy(vmu_tail, 1);
	sigma_tail = vec_copy(vsigma_tail, 1);
	mu_gamma = vec_copy(vmu_gamma, k);
	sigma_gamma = vec_copy(vsigma_gamma, k);
	mu_omega = vec_copy(vmu_omega, k);
	sigma_omega = vec_copy(vsigma_omega, k);
	mu_alpha = vec_copy(vmu_alpha, k);
	sigma_alpha = vec_copy(vsigma_alpha, k);
	mu_beta = vec_copy(vmu_beta, k);
	sigma_beta = vec_copy(vsigma_beta, k);
	mu_a = vec_copy(vmu_a, 1);
	sigma_a = vec_copy(vsigma_a, 1);
	mu_b = vec_copy(vmu_b, 1);
	sigma_b = vec_copy(vsigma_b, 1);
}


void memoryDeallocation(){
	_del(H); _del(Q); _del(R);
	_del(cholH);
	_del(cholH1);
	_del(MEs);
	_del(mH);
	_del(mMeanH);
	_del(y);
	vec_del(mu_tail);
	vec_del(sigma_tail);
	vec_del(mu_gamma);
	vec_del(sigma_gamma);
	vec_del(mu_omega);
	vec_del(sigma_omega);
	vec_del(mu_alpha);
	vec_del(sigma_alpha);
	vec_del(mu_a);
	vec_del(sigma_a);
	vec_del(mu_b);
	vec_del(sigma_b);
	vec_del(mu_beta);
	vec_del(sigma_beta);
}


double dssged(double *x, double *gamma, double delta, int k, int islog){
	int i;
	double y, mu, sigma, m1, aux, logDens;

	m1 = gammafn(2.0/delta) / sqrt( gammafn(1.0/delta) * gammafn(3.0/delta) );
	aux = pow( gammafn(3.0/delta)/gammafn(1.0/delta), delta/2.0 );

	// termos que nao dependem de i
	logDens = log(2.0) + log( gammafn(3.0/delta)/gammafn(1.0/delta) )/2 - log(2/delta) - log(gammafn(1.0/delta));
	logDens = ((double) k)*logDens;

	// termos que dependem de i
	for(i=0;i<k;i++){
		mu = (gamma[i] - 1/gamma[i]) * m1;
		sigma = sqrt( gamma[i]*gamma[i] + 1.0/(gamma[i]*gamma[i]) -1.0 - mu*mu );
		if(x[i] >= -mu/sigma)
			y = (x[i]*sigma + mu)/gamma[i];
		else
			y = (x[i]*sigma + mu)*gamma[i];

		logDens += log( gamma[i]*sigma/(1.0+ gamma[i]*gamma[i]) ) - aux*pow(fabs(y), delta);
	}

	if(islog)
		return logDens;
	else
		return exp(logDens);
}

void dssged_R(double *x, double *gamma, double *delta, int *k, int *islog, double *value){
	value[0] = dssged(x, gamma, delta[0], k[0], islog[0]);
}


double dsst(double *x, double *gamma, double v, int k, int islog){
	int i;
	double y, mu, sigma, m1, aux, logDens, dk;

	if(v>200) return(dssged(x, gamma, 2.0, k, islog));

	dk = (double) k;

	//m1 = 2*sqrt(v-2.0)*gammafn(0.5)*gammafn(v/2.0)/((v-1)*gammafn((v+1.0)/2.0));
	m1 = 0.5641896*sqrt(v-2.0)*gammafn((v-1)/2) / gammafn( v/2);

	// termos que nao dependem de i
	//logDens = dk*log(2.0) + log( gammafn((v+dk)/2.0)/gammafn(v/2.0) ) - dk*log(M_PI*(v-2.0))/2.0;
	logDens = dk*log(2.0) + lgammafn( (v+dk)/2.0 ) -lgammafn(v/2.0) - dk*log(M_PI*(v-2.0))/2.0;

	// termos que dependem de i
	aux=0.0;
	for(i=0;i<k;i++){
		mu = (gamma[i] - 1/gamma[i])*m1;
		sigma = sqrt( gamma[i]*gamma[i] + 1.0/(gamma[i]*gamma[i]) -1.0 - mu*mu );

		if(x[i] >= -mu/sigma)
			y = (x[i]*sigma + mu)/gamma[i];
		else
			y = (x[i]*sigma + mu)*gamma[i];

		logDens += log( gamma[i]*sigma/(1.0+ gamma[i]*gamma[i]) );
		aux += y*y;	// sum (x*)^2
	}

	logDens +=  -0.5*(v+dk)*log(1.0 + aux/(v-2.0));

	if(islog)
		return logDens;
	else
		return exp(logDens);
}

void dsst_R(double *x, double *gamma, double *v, int *k, int *islog, double *value){
	value[0] = dsst(x, gamma, v[0], k[0], islog[0]);
}


double logLikelihood(double *omega, double *alpha, double *beta, double a, double b, double *gamma, double tail){
	int cholTest;
	int i,j;
	int t;
	double ht1[k]; // h_{ii,t-1}, i=1,...,k
	double hiit;
	double z[k]; // z_t^* vector
	double value=0.0;
	double Inf = 999999999999999999;
	double logDetH;

	mat_zero(H, k);
	mat_zero(cholH, k);
	mat_zero(cholH1, k);
	mat_zero(Q, k);
	mat_zero(R, k);

	// compute the R matrix
	for(i = 0; i < k; i++){
		//hiit = ( y[0][i]*y[0][i] + y[1][i]*y[1][i] + y[2][i]*y[2][i] + y[3][i]*y[3][i] + y[4][i]*y[4][i] )/5 ; // omega[i];// /(1.0-beta[i]); // H_{ii,1} // modificado na versao 2.1
		hiit = omega[i] / (1.0 - alpha[i] - beta[i]); // modificado na versao 2.3 (06/07/2021)
		MEs[i][0]  = y[0][i]/sqrt(hiit);
		for(t = 1; t < n; t++){
			hiit = omega[i] + alpha[i]*y[t-1][i]*y[t-1][i] + beta[i]*hiit; // H_{ii,t}
			MEs[i][t]  = y[t][i]/sqrt(hiit); // Standard Errors
		}
	}
	mcov(n, k, MEs, R); // compute the R matrix

	//*************** Initialization ****************//
	//  H_{1} and Q_{1}
	for(i = 0; i < k; i++){
		Q[i][i] = R[i][i]; //(1.0-a-b)*R[i][i];// /(1.0-b);
		//ht1[i] =  ( y[0][i]*y[0][i] + y[1][i]*y[1][i] + y[2][i]*y[2][i] + y[3][i]*y[3][i] + y[4][i]*y[4][i] )/5; //omega[i];// /(1.0-beta[i]); // modificado na versao 2.1
		ht1[i] = omega[i] / (1.0 - alpha[i] - beta[i]); // modificado na versao 2.3 (06/07/2021)
 		H[i][i] = ht1[i];
		for(j=i+1; j<k; j++){
			Q[i][j] = R[i][j];// (1.0-a-b)*R[i][j];// /(1.0-b); // modificado na versao 2.1
			Q[j][i] = Q[i][j];
			H[i][j] = R[i][j];     // Note that it is not used
			H[j][i] = H[i][j];     // Note that it is not used
		}
	}

	matrix_to_vector(k, k, H, mH[0]); // save the H_0 matrix
	//***************************************************//

	// ****************** main looping ******************//
	t=1; // This is equivalent to t=2 in the model
	while(t<n){
		for(i=0;i<k;i++){
			Q[i][i] = (1-a-b)*R[i][i]  +  a*y[t-1][i]*y[t-1][i]/H[i][i]  +  b*Q[i][i]; //Q_{ii,t}
			H[i][i] = omega[i] + alpha[i]*y[t-1][i]*y[t-1][i] + beta[i]*H[i][i]; //H_{ii,t}
		}
		for(i=0;i<k;i++){
			for(j=i+1; j<k; j++){
				Q[i][j] = (1-a-b)*R[i][j] + a*y[t-1][i]*y[t-1][j]/sqrt( ht1[i]*ht1[j] ) + b*Q[i][j]; //Q_{ij,t} (Note that it is function of H_{ii,t-1} and H_{jj,t-1})
				Q[j][i] = Q[i][j]; //Q_{ji,t} (Note that it is not used.)
				H[i][j] = Q[i][j] * sqrt( (H[i][i]*H[j][j]) / (Q[i][i]*Q[j][j]) ); //H_{ij,t}
				H[j][i] = H[i][j]; //H_{ji,t}
			}
			ht1[i] = H[i][i]; //H_{ii,t} (Note that H_{ii,t-1} will not be used forward in this loop)
		}

		matrix_to_vector(k, k, H, mH[t]); // save the H_t matrix


		cholTest = choldc2(H,  cholH, k);
		logDetH = 2 * mat_triangular_det(cholH, k, 1); // log |H_t| = 2*log( |H_t^{1/2}| )

		if( ( logDetH < -100 ) || cholTest == 0 ){
			//Rprintf("\nH is a singular matrix or H is not positive definite.\n logDetH = %.10lf\n", logDetH);
			//RprintMatrix(k, k, H);
			return -Inf;
		}

		mat_inv_L(cholH, cholH1, k);

		mult_mat_vector(cholH1, y[t], z, k, k); // z_t = H_{t}^{-1/2} y_t

		value += -0.5 * logDetH; //  -1/2 log(|H_t|)

		if(errorDist==1) //ssnorm
			value += dssged(z, gamma, 2.0, k, 1); // + log p_e ( H_{t}^{-1/2} y_t )
		if(errorDist==2) //sst
			value += dsst(z, gamma, tail, k, 1); // + log p_e ( H_{t}^{-1/2} y_t )
		if(errorDist==3) //ssged
			value += dssged(z, gamma, tail, k, 1); // + log p_e ( H_{t}^{-1/2} y_t )

		t++;
	}

	logLikelihood_last = value; //used for DIC

	return(value);
}

void logLikelihood_R(double *vY, double *omega, double *alpha, double *beta,
						double *a, double *b, double *gamma, double *tail, int *verrorDist,
						int *vn, int *vk, double *value, double *vmH){
	n = vn[0];
	k = vk[0];
	errorDist = verrorDist[0];

	H=mat_new(k);  Q=mat_new(k); R=mat_new(k);
	cholH=mat_new(k);
	cholH1=mat_new(k); //Cholesky decomposition of H1
	matrix_new(k, n, &MEs); //allocation of errors matrix (Note that the order is k X n )

	matrix_new(n, k, &y); // allocation of the data matrix
	vector_to_matriz(n, k, vY, y);

	matrix_new(n, k*k, &mH); // matrix to save the H_t's for t=1,...,n

	value[0] = logLikelihood(omega, alpha, beta, a[0], b[0], gamma, tail[0]);

	matrix_to_vector(n, k*k, mH, vmH);

	_del(H); _del(Q); _del(R); _del(cholH);
	_del(cholH1);
	_del(MEs);
	_del(mH);
	_del(y);
}


/*
	Input: one markov chain for parameters
	Output: loglikelihood mean (also update the mean H)
	Arguments:
		mcmc: matrix of parameters (original scale)
		n_mcmc: number of lines of mcmc
*/
double loglike_matrix(double **mcmc, int n_mcmc){
	int i,j;
	double meanLogLike, *gamma, *omega, *alpha, *beta, a, b, tail;

	gamma=vec_new(n_mcmc); omega=vec_new(n_mcmc); alpha=vec_new(n_mcmc); beta=vec_new(n_mcmc);

	matrix_zero(mMeanH, n, k*k);
	matrix_zero(mH, n, k*k);

	meanLogLike = 0.0;
	for(j=0; j<n_mcmc; j++){

		tail = mcmc[j][0];

		for(i=1; i<=k; i++){
			gamma[i-1] = mcmc[j][4*(i-1)+1];
			omega[i-1] = mcmc[j][4*(i-1)+2];
			alpha[i-1] = mcmc[j][4*(i-1)+3];
			beta[i-1] = mcmc[j][4*i];
		}

		a = mcmc[j][4*k +1];
		b = mcmc[j][4*k +2];

		meanLogLike += logLikelihood(omega, alpha, beta, a, b, gamma, tail) / ((double) n_mcmc );
		matrix_sum(mMeanH, mMeanH, mH, n, k*k);
	}

	matrix_mult_cte(mMeanH, mMeanH, 1.0/((double) n_mcmc), n, k*k);

	vec_del(gamma); vec_del(omega); vec_del(alpha); vec_del(beta);

	return( meanLogLike );
}

void loglike_matrix_R(double *vY, int *vn, int *vk, double *vmcmc, int *n_mcmc, int *verrorDist, double *meanLogLike, double *vmMeanH){
	double **mcmc;
	n = vn[0];
	k = vk[0];
	errorDist = verrorDist[0];

	H=mat_new(k);  Q=mat_new(k); R=mat_new(k);
	cholH=mat_new(k);
	cholH1=mat_new(k); //Cholesky decomposition of H1
	matrix_new(k, n, &MEs); //allocation of errors matrix (Note that the order is k X n )

	matrix_new(n, k, &y); // allocation of the data matrix
	vector_to_matriz(n, k, vY, y);

	matrix_new(n_mcmc[0], (4*k) +3, &mcmc);
	vector_to_matriz(n_mcmc[0], (4*k) +3, vmcmc, mcmc);

	matrix_new(n, k*k, &mH); // matrix to save the H_t's for t=1,...,n
	matrix_new(n, k*k, &mMeanH); // matrix to save the mean of H_t's for t=1,...,n

	meanLogLike[0] = loglike_matrix(mcmc, n_mcmc[0]);

	//RprintMatrix(n, k*k, mMeanH);

	matrix_to_vector(n, k*k, mMeanH, vmMeanH);

	_del(mcmc);
	_del(H); _del(Q); _del(R); _del(cholH);
	_del(cholH1);
	_del(MEs);
	_del(mH); _del(mMeanH);
	_del(y);
}


double logPosterior(double *omega, double *alpha, double *beta, double a, double b, double *gamma, double tail){
	int i;
	double logLik, logPrior=0.0;

	double Inf = 999999999999999999;

	if(a+b>0.989) return -Inf;
	for(i = 0; i < k; i++)
		if(alpha[i]+beta[i]>0.989) return -Inf;


	logLik = logLikelihood(omega, alpha, beta, a, b, gamma, tail);

	if(k>1)
		logPrior = dnorm(a, mu_a[0], sigma_a[0], 1) + dnorm(b, mu_b[0], sigma_b[0], 1);
	else
		logPrior = 0.0;

	if(errorDist > 1)
		logPrior += dnorm(tail, mu_tail[0], sigma_tail[0], 1);

	for(i=0;i<k;i++)
		logPrior +=  dnorm(omega[i], mu_omega[i], sigma_omega[i], 1) + dnorm(alpha[i], mu_alpha[i], sigma_alpha[i], 1)
				   + dnorm(beta[i], mu_beta[i], sigma_beta[i], 1) + dnorm(gamma[i], mu_gamma[i], sigma_gamma[i], 1);

	//Rprintf("\nlogLik: %f, logPrior: %f\n", logLik, logPrior);

	return(logLik + logPrior);
}


void logPosterior_R(double *vY, double *omega, double *alpha, double *beta,
						double *a, double *b, double *gamma, double *tail, int *verrorDist,
						int *vn, int *vk, double *vmH,
						double *vmu_omega,
						double *vsigma_omega,
						double *vmu_alpha,
						double *vsigma_alpha,
						double *vmu_beta,
						double *vsigma_beta,
						double *vmu_a,
						double *vsigma_a,
						double *vmu_b,
						double *vsigma_b,
						double *vmu_gamma,
						double *vsigma_gamma,
						double *vmu_tail,
						double *vsigma_tail,
						double *value){
	n=vn[0];
	k=vk[0];
	errorDist=verrorDist[0];

	H=mat_new(k); Q=mat_new(k); R=mat_new(k); //P=mat_new(k[0]); L=mat_new(k[0]); U=mat_new(k[0]);
	cholH = mat_new(k);
	cholH1=mat_new(k); //Cholesky decomposition of H1
	matrix_new(k, n, &MEs); //allocation of errors matrix (Note that the order is k X n )

	matrix_new(n, k, &y); // allocation of the data matrix
	vector_to_matriz(n, k, vY, y);

	matrix_new(n,k*k, &mH); // matrix to save the H_t's for t=1,...,n

	value[0] = logLikelihood(omega, alpha, beta, a[0], b[0], gamma, tail[0]);

	matrix_to_vector(n, k*k, mH, vmH);

	mu_omega = vmu_omega;
	sigma_omega = vsigma_omega;
	mu_alpha = vmu_alpha;
	sigma_alpha = vsigma_alpha;
	mu_beta = vmu_beta;
	sigma_beta = vsigma_beta;
	mu_a = vmu_a;
	sigma_a = vsigma_a;
	mu_b = vmu_b;
	sigma_b = vsigma_b;
	mu_gamma = vmu_gamma;
	sigma_gamma = vsigma_gamma;
	mu_tail = vmu_tail;
	sigma_tail = vsigma_tail;

	value[0] = logPosterior(omega, alpha, beta, a[0], b[0], gamma, tail[0]);

	_del(H); _del(Q); _del(R);
	_del(cholH);
	_del(cholH1);
	_del(MEs);
	_del(mH);
	_del(y);
}



/* 	Real line transformation for all parameters:

	phi[0]    			=  log( tail - 2 )   if errorDist=2          // t-student distribution
						=  log( tail )       if errorDist=3          // ged distribution (or normal distribution)

	phi[1:4]  			=  ( log(gamma_1), log(omega_1), log(alpha_1/(1-alpha_1)), log(beta_1/(1-beta_1)) )'
	phi[5:8]  			=  ( log(gamma_2), log(omega_2), log(alpha_2/(1-alpha_2)), log(beta_2/(1-beta_2)) )'
	...
	phi[(4(k-1)+1):4k]  =  ( log(gamma_k), log(omega_k), log(alpha_k/(1-alpha_k)), log(beta_k/(1-beta_k)) )'

	phi[4k+1]			=  log( a/(1-a) )
	phi[4k+2]			=  log( b/(1-b) )
*/
void real_to_original_scale(double *phi, double *omega, double *alpha, double *beta,
							double *a, double *b, double *gamma, double *tail, int *k, int *errorDist){
	int i;

	if(*errorDist == 1 || *errorDist == 3)
		tail[0] = exp(phi[0]);
	if(*errorDist == 2)
		tail[0] = exp(phi[0])+2.0;

	for(i=1; i<=k[0]; i++){
		gamma[i-1] = exp( phi[ 4*(i-1)+1 ] );
		omega[i-1] = exp( phi[ 4*(i-1)+2 ] );
		alpha[i-1] = exp( phi[ 4*(i-1)+3 ] ) / ( 1.0+exp( phi[ 4*(i-1)+3 ] ) );
		beta[i-1]  = exp( phi[ 4*i ] ) / ( 1.0+exp( phi[ 4*i ] ) );
	}
	a[0]  = exp( phi[ 4*k[0] +1 ] ) / ( 1.0+exp( phi[ 4*k[0] +1 ] ) );
	b[0]  = exp( phi[ 4*k[0] +2 ] ) / ( 1.0+exp( phi[ 4*k[0] +2 ] ) );
}

/* 	Real line transformation for all parameters:

	phi[0]    			=  log( tail - 2 )   if errorDist=2          // t-student distribution
						=  log( tail )       if errorDist=3          // ged distribution (or normal distribution)

	phi[1:4]  			=  ( log(gamma_1), log(omega_1), log(alpha_1/(1-alpha_1)), log(beta_1/(1-beta_1)) )'
	phi[5:8]  			=  ( log(gamma_2), log(omega_2), log(alpha_2/(1-alpha_2)), log(beta_2/(1-beta_2)) )'
	...
	phi[(4(k-1)+1):4k]  =  ( log(gamma_k), log(omega_k), log(alpha_k/(1-alpha_k)), log(beta_k/(1-beta_k)) )'

	phi[4k+1]			=  log( a/(1-a) )
	phi[4k+2]			=  log( b/(1-b) )
*/
void original_to_real_scale(double *phi, double *omega, double *alpha, double *beta,
							double *a, double *b, double *gamma, double *tail, int *k, int *errorDist){

	int i;

	if(*errorDist == 1 || *errorDist == 3)
		phi[0] = log(tail[0]);
	if(*errorDist == 2)
		phi[0] = log(tail[0]-2.0);

	for(i=1; i<=k[0]; i++){
		phi[ 4*(i-1)+1 ] = log(gamma[i-1]);
		phi[ 4*(i-1)+2 ] = log(omega[i-1]);
		phi[ 4*(i-1)+3 ] = log(alpha[i-1]/(1.0-alpha[i-1]));
		phi[ 4*i ] = log(beta[i-1]/(1.0-beta[i-1]));
	}

	phi[ 4*k[0] +1 ] = log(a[0]/(1.0-a[0]));
	phi[ 4*k[0] +2 ] = log(b[0]/(1.0-b[0]));
}



double logPosterior_phi(double *phi){
	int i;
	double logPost, logJacob;
	double omega[k], alpha[k], beta[k], a, b, gamma[k], tail=3;

	real_to_original_scale(phi, omega, alpha, beta,
							&a, &b,  gamma,  &tail, &k, &errorDist);

	logPost = logPosterior(omega, alpha, beta, a, b, gamma, tail);


	if(errorDist == 1) logJacob = 0.0;
	else 	logJacob = phi[0];


	for(i=1; i <= k; i++){
		logJacob += phi[4*(i-1)+1];  // log d gamma_i )/ d phi
		logJacob += phi[4*(i-1)+2];  // log d omega_i / d phi
		logJacob += phi[4*(i-1)+3] - 2*log( 1.0 + exp( phi[4*(i-1)+3] ) ); // log d alpha_i / d phi
		logJacob += phi[ 4*i ] - 2*log( 1.0 + exp(phi[4*i]) );  // log d beta_i / d phi
	}


	if( k > 1 )
		logJacob += phi[4*k+1] - 2*log( 1.0 + exp(phi[4*k+1]) ) + phi[4*k+2] - 2*log( 1.0 + exp(phi[4*k+2]) ); // log d a /d phi + // log d b /d phi

	return( logPost+logJacob );
}

void logPosterior_phi_R(double *phi, double *value){
	value[0] = logPosterior_phi(phi);
}


void MH_oneDimension(double *phi,
					double *sd_phi_sim, // vector of standard deviations for generate candidates in the simulation
					int *n_sim, // number of required simulations
					double *vMC // monte carlo simulation
				)
{
	int s,i,start;
	int i_vMC; //index of vMC
	int n_par; // number of parameters
	double logPoster, logPoster_new;
	double phi_ti; // phi_{t,i}, i=1,...,n_par
	double uu, vv;
	double **mH_actual;

	if(k == 1) n_par=5; // five
	else n_par = 4*(k) + 3; // 4k + 3

	// ****************** memory allocation *********************** //
	matrix_new(n,k*k, &mH_actual);
	// ************************************************************* //

	logPoster = logPosterior_phi(phi);
	logPoster_new = logPoster;
	logLikelihood_old = logLikelihood_last;
	logLikelihood_mean = 0.0 ; // used for DIC information
	matrix_zero(mMeanH, n, k*k);
	matrix_copy(mH_actual, mH, n, k*k);

	start=0;
	i_vMC=0;
	s=0;
	while(s < *n_sim){

		if((s+1)%100 == 0 && print==1) Rprintf("Simulation number = %d\n", s+1);

		if(errorDist==1){ start=1; vMC[i_vMC]=log(2.0); i_vMC++; }

		for(i=start; i< n_par; i++){

			GetRNGstate();
			vv = norm_rand();
			PutRNGstate();

			phi_ti = phi[i];
			phi[i] = phi_ti + sd_phi_sim[i]*vv;

			logPoster_new = logPosterior_phi(phi);

			GetRNGstate();
			uu = unif_rand();
			PutRNGstate();

			if( uu < exp( logPoster_new -logPoster ) ){
				logPoster = logPoster_new;
				matrix_copy(mH_actual, mH, n, k*k);
				logLikelihood_old = logLikelihood_last;
			}else{
				phi[i] = phi_ti;
			}

			vMC[i_vMC] = phi[i];
			i_vMC++;
		}

		matrix_sum(mMeanH, mMeanH, mH_actual, n, k*k);
		logLikelihood_mean += logLikelihood_old/((double) n_sim[0]) ;

		s++;
	}

	matrix_mult_cte(mMeanH, mMeanH, 1.0/((double) n_sim[0]), n, k*k);

	//RprintMatrix(n, k*k, mMeanH);

	// ****************** memory deallocation *************** //
	_del(mH_actual);
	// ***************************************************** //
}

void rMultNorm(double *mean, mat chol_cov, double *out, _dim){
	double z[n];

	for_i{	GetRNGstate(); z[i] = norm_rand(); PutRNGstate(); }

	mult_mat_vector( chol_cov, z, out, n, n);

	for_i out[i] += mean[i];

	return;
}

/*
void rMultNorm(double *mu, double **sigma, double *x, int l){

	 int cont;
	 double u1,u2,m1,m2,aux[l], aux2[l];

	 cont = 0;
	 while(cont < l){ // gera o vetor aleatorio da normal padrao na variavel "aux"
		 u1 = ( (double)rand() + 0.0001 )/( (double)RAND_MAX + 0.0002);
		 u2 = ( (double)rand() + 0.0001 )/( (double)RAND_MAX + 0.0002);
		 m1 = -2.0*log(u1);
		 m1 = sqrt(m1);
		 m2 = 2.0*3.1415926535897932*u2;
		 aux[cont++] = m1*cos(m2);
		 if(cont < l)
			aux[cont++] = m1*sin(m2);
	 }

	 mult_mat_vector( sigma, aux, aux2, l, l);

	 for(cont = 0; cont<l; cont++)
		x[cont] = aux2[cont] + mu[cont];

}*/


void MH_oneBlock( double *phi,
				  double *vChol_Cov_phi_sim, // Cholesky decomposition of the covariance matrix for generate candidates
				  int *n_sim, // number of required simulations
				  double *vMC // monte carlo simulation
				)
{
	int s,i;
	int i_vMC; //index of vMC
	int n_par; // number of parameters
	double logPoster, logPoster_new;
	double uu;

	if(k == 1) n_par=5; // five
	else n_par = 4*(k) + 3; // 4k + 3

	vec phi_new;
	mat  Chol_Cov_phi_sim;
	double **mH_actual;

	// ****************** memory allocation *********************** //
	phi_new = vec_copy(phi, 4*(k) + 3);
	Chol_Cov_phi_sim = mat_new(n_par);
	vector_to_matriz(n_par, n_par, vChol_Cov_phi_sim, Chol_Cov_phi_sim);
	matrix_new(n, k*k, &mH_actual);

	//Rprintf("\nMatrixCholCovSim: ");
	//RprintMatrix(n_par, n_par, Chol_Cov_phi_sim);
	// ************************************************************* //

	logPoster = logPosterior_phi(phi);
	logPoster_new = logPoster;
	logLikelihood_old = logLikelihood_last;
	logLikelihood_mean = 0.0; // used for DIC information
	matrix_zero(mMeanH, n, k*k);
	matrix_copy(mH_actual, mH, n, k*k);

	i_vMC=0;
	s=0;
	while(s < *n_sim){

		if((s+1)%100 == 0 && print==1) Rprintf("Simulation number = %d\n", s+1);

		rMultNorm(phi, Chol_Cov_phi_sim, phi_new, n_par);

		if(errorDist==1){ phi_new[0] = log(2.0); }

		logPoster_new = logPosterior_phi(phi_new);

		GetRNGstate();
		uu = unif_rand();
		PutRNGstate();

		if( uu < exp( logPoster_new -logPoster ) ){
			logPoster = logPoster_new;

			matrix_copy(mH_actual, mH, n, k*k);
			logLikelihood_old = logLikelihood_last;

			for(i=0; i<n_par; i++) phi[i] = phi_new[i];
		}

		for(i=0; i<n_par; i++){
			vMC[i_vMC] = phi[i];
			i_vMC++;
		}

		matrix_sum(mMeanH, mMeanH, mH_actual, n, k*k);
		logLikelihood_mean += logLikelihood_old/((double) n_sim[0]) ;

		s++;
	}

	matrix_mult_cte(mMeanH, mMeanH, 1.0/((double) n_sim[0]), n, k*k);


	// ****************** memory deallocation *************** //
	vec_del(phi_new);
	_del(Chol_Cov_phi_sim);
	_del(mH_actual);
	// ***************************************************** //
}
