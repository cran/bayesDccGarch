
#ifndef _MAT_OPE_H_
#define _MAT_OPE_H_


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>


#define foreach(a, b, c) for (int a = b; a < c; a++)
#define for_i foreach(i, 0, n)
#define for_j foreach(j, 0, n)
#define for_k foreach(k, 0, n)
#define for_ij for_i for_j
#define for_ijk for_ij for_k
#define dim_ int n
#define swap_(x, y) { typeof(x) tmp = x; x = y; y = tmp; }
#define sum_k_(a, b, c, s) { s = 0; foreach(k, a, b) s+= c; }


typedef double **mat;
typedef double *vec;


void vec_zero(vec x, dim_){
	for_i x[i] = 0.0;
}

vec vec_new(dim_)
{
	vec x = malloc(sizeof(double) * n);

	vec_zero(x, n);

	return x;
}

void vec_del(vec x) { free(x); }

vec vec_copy(vec v, dim_)
{
	vec x = vec_new(n);
	for_i x[i] = v[i];
	return x;
}


void matrix_new(int nRows, int nColluns, double ***M){
    int i;//,j;

	// Allocation
    (*M) = (double**) malloc(nRows * sizeof (double*));
    for (i = 0; i < nRows; i++){
        (*M)[i] = (double*) malloc(nColluns * sizeof (double));
    }

	//Initialization
	/*
	for (i = 0; i < nRows; i++)
		for (j = 0; j < nColluns; j++)
			*M[i][j] = 0.0;
    */
	return;
}

void matrix_zero(double **x, int n, int m){
	int i,j;
	for(i=0;i<n;i++){ for(j=0;j<m;j++){ x[i][j] = 0.0; } }
}

// a = b+c
void matrix_sum(double **a, double **b, double **c, int n, int m){
	int i,j;
	for(i=0;i<n;i++){
	  for(j=0;j<m;j++){
	    a[i][j] = b[i][j] + c[i][j];
	  }
	}
}

// a = b
void matrix_copy(double **a, double **b, int n, int m){
	int i,j;
	for(i=0;i<n;i++){
	  for(j=0;j<m;j++){
	    a[i][j] = b[i][j];
	  }
	}
}

// a = cte*b
void matrix_mult_cte(double **a, double **b, double cte, int n, int m){
	int i,j;
	for(i=0;i<n;i++){
	  for(j=0;j<m;j++){
	    a[i][j] = cte*b[i][j];
	  }
	}
}


#define _zero(a) mat_zero(a, n)
void mat_zero(mat x, int n) { for_ij x[i][j] = 0; }

#define _new(a) a = mat_new(n)
mat mat_new(dim_)
{
	mat x = malloc(sizeof(double*) * n);
	x[0]  = malloc(sizeof(double) * n * n);

	for_i x[i] = x[0] + n * i;
	_zero(x);

	return x;
}

#define _copy(a) mat_copy(a, n)
mat mat_copy(void *s, dim_)
{
	mat x = mat_new(n);
	for_ij x[i][j] = ((double (*)[n])s)[i][j];
	return x;
}

#define _trans(a,b) mat_trans(a,b, n)
void mat_trans(mat A, mat At, dim_)
{
	for_ij At[i][j] = A[j][i];
}

#define _del(x) mat_del(x)
void mat_del(mat x) { free(x[0]); free(x); }


#define _mul(a, b) mat_mul(a, b, n)
mat mat_mul(mat a, mat b, dim_)
{
	mat c = _new(c);
	for_ijk c[i][j] += a[i][k] * b[k][j];
	return c;
}

#define _mul2(a, b, c) mat_mul2(a, b, c, n)
void mat_mul2(mat a, mat b, mat c, dim_)
{
	_zero(c);
	for_ijk c[i][j] += a[i][k] * b[k][j];
}

/*  M * V = MxV, where M(nxq), V(qx1) and  MxV(nx1)*/
void mult_mat_vector(mat M, double *V, double *MxV, int n, int q){
	int i,j;
	double aux;
	for(i=0;i<n;i++){
		aux=0.0;
		for(j=0;j<q;j++)
			aux += M[i][j]*V[j];
		MxV[i] = aux;
	}
}


#define _pivot(a, b) mat_pivot(a, b, n)
void mat_pivot(mat a, mat p, dim_)
{
	for_ij { p[i][j] = (i == j); }
	for_i  {
		int max_j = i;
		foreach(j, i, n)
			if (fabs(a[j][i]) > fabs(a[max_j][i])) max_j = j;

		if (max_j != i)
			for_k { swap_(p[i][k], p[max_j][k]); }
	}
}

//convert vector to matrix by row
void vector_to_matriz(int n, int q, double *V, double **M){
	int i,j;
	for (i=0; i<n; i++){
		for (j=0; j<q; j++){
			M[i][j] = V[i*q+j];
		}
	}
}

//convert matrix to vector by row
void matrix_to_vector(int n, int q, double **M, double *V){
	int a, b;
	for (a=0; a<n; a++){
		for (b=0; b<q; b++){
			V[a*q + b] = M[a][b];
		}
	}
}



#define _PLU(a, p, l, u) mat_PLU(a, p, l, u, n)
void mat_PLU(mat A, mat P, mat L, mat U, dim_)
{
	_zero(L); _zero(U);
	_pivot(A, P);

	mat Aprime = _mul(P, A);

	for_i  { L[i][i] = 1; }
	for_ij {
		double s;
		if (j <= i) {
			sum_k_(0, j, L[j][k] * U[k][i], s)
			U[j][i] = Aprime[j][i] - s;
		}
		if (j >= i) {
			sum_k_(0, i, L[j][k] * U[k][i], s);
			L[j][i] = (Aprime[j][i] - s) / U[i][i];
		}
	}

	_del(Aprime);
}


void mat_PLU_v(double *vA, double *vP, double *vL, double *vU, int *n)
{
	mat A, L, U, P;
	A = mat_new(*n);
	L = mat_new(*n);
	U = mat_new(*n);
	P = mat_new(*n);
	vector_to_matriz(*n, *n, vA, A);
	vector_to_matriz(*n, *n, vL, L);
	vector_to_matriz(*n, *n, vU, U);
	vector_to_matriz(*n, *n, vP, P);

	mat_PLU(A, P, L, U, *n);

	matrix_to_vector(*n, *n, A, vA);
	matrix_to_vector(*n, *n, L, vL);
	matrix_to_vector(*n, *n, U, vU);
	matrix_to_vector(*n, *n, P, vP);

	mat_del(A); mat_del(L); mat_del(U); mat_del(P);
}


void mat_solve_Lzb(mat L, double *z, double *b, dim_){
	int i,j;
	z[0] = b[0]/L[0][0];
	for(i=1; i<n; i++){
		z[i] = b[i];
		for(j=0; j<i; j++)
			z[i] -= L[i][j]*z[j];
		z[i] /= L[i][i];
	}
}

void mat_solve_Uzb(mat U, double *z, double *b, dim_){
	int i,j;
	z[n-1] = b[n-1]/U[n-1][n-1];
	for(i=n-2; i>=0; i--){
		z[i] = b[i];
		for(j=i+1; j<n; j++)
			z[i] -= U[i][j]*z[j];
		z[i] /= U[i][i];
	}
}

// // Compute the Inverse of a matrix AnXn using LU decomposition
#define _inv(b, p, l, u) mat_inv(b, p, l, u, n)
void mat_inv(mat Ainv, mat P, mat L, mat U, dim_){

	if(n==1){ Ainv[0][0] = 1.0/U[0][0]; return;}

	_zero(Ainv);

	mat Aaux,Pt;
	_new(Aaux); _new(Pt);

	double z[n], b[n], I[n];

	for_i{

		for_j{ z[j]=0.0; b[j]=0.0; I[j]=0.0; }

		I[i] = 1.0;

		mat_solve_Lzb(L, z, I, n);
		mat_solve_Uzb(U, b, z, n);

		for_j{ Aaux[j][i] = b[j]; }
	}

	_trans(P, Pt);
	_mul2(Aaux, Pt, Ainv);

	_del(Aaux); _del(Pt);

}


#define _inv2(a) mat_inv2(a, n)
mat mat_inv2(mat A, dim_){

	mat P, L, U, Ainv;

	_new(P); _new(L); _new(U); _new(Ainv);

	_PLU(A, P, L, U);

	_inv(Ainv, P, L, U);

	_del(P); _del(L);  _del(U);

	return(Ainv);
}


void mat_inv_v(double *vA, double *vAinv, int *n){
	mat A, Ainv;

	A = mat_new(*n);

	vector_to_matriz(*n, *n, vA, A);

	Ainv = mat_inv2( A, *n);

	matrix_to_vector(*n, *n, Ainv, vAinv);

	mat_del(A);  mat_del(Ainv);
}



// inverse of a lower triangular matrix (L)
void mat_inv_L(mat L, mat L1, dim_){
	int i,j,k;
	double aux;

	for(i=0; i<n; i++){
		L1[i][i] = 1.0/L[i][i];
		for(j=0; j<i; j++){
			aux = 0.0;
			for(k=j; k<i; k++){
				aux -= L[i][k]*L1[k][j];
			}
			aux /= L[i][i];
			L1[i][j] = aux;
			L1[j][i] = 0.0;
		}
	}

}

void mat_inv_L_v(double *vL, double *vLinv, int *n){
	mat L, L1;

	L = mat_new(*n);
	L1 = mat_new(*n);

	vector_to_matriz(*n, *n, vL, L);
	vector_to_matriz(*n, *n, vLinv, L1);

	mat_inv_L(L, L1, *n);

	matrix_to_vector(*n, *n, L1, vLinv);

	_del(L);  _del(L1);
}

double mat_triangular_det(mat A, dim_, int islog){
	double logDet;

	logDet=0.0;
	for_i logDet += log(A[i][i]);

	if(islog)
		return(logDet);
	else
		return(exp(logDet));
}

void mat_triangular_det_v(double *vL, double *detL, int *n){
	mat L;

	L = mat_new(*n);

	vector_to_matriz(*n, *n, vL, L);

	*detL = mat_triangular_det(L, *n, 0);

	_del(L);
}

// // Compute the Determinant of a matrix AnXn using LU decomposition
#define _det(p, l, u) mat_det(p, l, u, n)
double mat_det(mat P, mat L, mat U, dim_){

	if(n==1) return(U[0][0]);

	double s,mult;

	mult=1.0;
	s=0;
	for_i{
		if(P[i][i]== 0.0){s += 1.0;}
		mult *= U[i][i];
	}

	mult *= pow(-1, s/2.0);

	return(mult);
}


#define _det2(a) mat_det2(a, n)
double mat_det2(mat A, dim_){
	double det;
	mat P, L, U;

	_new(P); _new(L); _new(U);

	_PLU(A, P, L, U);

	det = _det(P, L, U);

	_del(P); _del(L);  _del(U);

	return(det);
}



/*Given a positive-definite symmetric matrix a[0..n-1][0..n-1], this routine constructs its Cholesky
decomposition, A = L · LT . On input, only the upper triangle of a need be given; it is not
modified. The Cholesky factor L is returned in the lower triangle of L.*/
int choldc(mat a,  mat L, dim_)
{
	if(n==1){ L[0][0] = sqrt(a[0][0]); return 1; }
	//void nrerror(char error_text[]);
	int i,j,k;
	double sum;

	for (i=1;i<=n;i++) {
		for (j=i;j<=n;j++) {
			for (sum=a[i-1][j-1],k=i-1;k>=1;k--) sum -= a[i-1][k-1]*a[j-1][k-1];
			if (i == j) {
				if (sum <= 0.0){ //a, with rounding errors, is not positive definite.
					Rprintf("\n\ncholdc failed!\n\n");
					return 0;
				}
				L[i-1][i-1]=sqrt(sum);
			} else a[j-1][i-1]=sum/L[i-1][i-1];
		}
	}

	for(i=0;i<n;i++){
		for (j=0;j<i;j++){
			L[i][j] = a[i][j];
			L[j][i] = 0.0;
			a[i][j] = a[j][i];
		}
	}
	return 1;
}



/*Given a positive-definite symmetric matrix a[0..n-1][0..n-1], this routine constructs its Cholesky
decomposition, A = L · LT . On input, only the upper triangle of a need be given; it is not
modified. The Cholesky factor L is returned in the lower triangle of L.*/
int choldc2(mat a,  mat L, dim_)
{
	if(n==1){ L[0][0] = sqrt(a[0][0]); return 1; }
	//void nrerror(char error_text[]);
	int i,j,k;
	double sum;

	for (i=1;i<=n;i++) {
		for (j=i;j<=n;j++) {

			for (sum=a[i-1][j-1],k=i-1;k>=1;k--) sum -= L[i-1][k-1]*L[j-1][k-1];

			if (i == j) {
				if (sum <= 0.0){ //a, with rounding errors, is not positive definite.
					Rprintf("\n\ncholdc failed!\n\n");
					return 0;
				}
				L[i-1][i-1] = sqrt(sum);
			}else{
				L[j-1][i-1] = sum/L[i-1][i-1];
				L[i-1][j-1] = 0.0;
			}

		}
	}

	return 1;
}

void choldc2_v(double *va,  double *vL, int *n, int *test){
	mat L, a;

	L = mat_new(*n);
	a = mat_new(*n);

	vector_to_matriz(*n, *n, va, a);
	vector_to_matriz(*n, *n, vL, L);

	*test = choldc2(a,  L, *n);

	matrix_to_vector(*n, *n, L, vL);

	_del(L);  _del(a);
}


// Compute the covariance between the vectors X and Y
double cov(int n, double *X, double *Y){
	int i;
	double media_XY, media_X, media_Y;

	media_XY = 0.0; media_X = 0.0; media_Y = 0.0;
	for(i=0;i<n;i++){
		media_XY += X[i]*Y[i];
		media_X += X[i];
		media_Y += Y[i];
	}

	media_XY = media_XY / ( (double) n);
	media_X = media_X / ( (double) n);
	media_Y = media_Y / ( (double) n);

	return( media_XY - media_X*media_Y  );
}

// Compute the covariance matrix for the rows of the matrix X
void mcov(int n, int p, double **X, mat result){
	int i,j;

	for(i=0;i<p;i++)
		for(j=i;j<p;j++){
			result[i][j] = cov(n, X[i], X[j]);
			result[j][i] = result[i][j];
		}

	return;
}


void RprintVector(int n, double *V){
	int i;
	for(i = 0; i < n; i++){
			Rprintf("%.5f, ", V[i]);
		Rprintf("\n");
	}
}

void RprintMatrix(int n, int p, double **M){
	int i,j;
	for(i = 0; i < n; i++){
		for(j = 0; j < p; j++)
			Rprintf("%.5f, ", M[i][j]);
		Rprintf("\n");
	}
}

#endif
