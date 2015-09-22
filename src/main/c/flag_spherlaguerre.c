 // FLAG package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "flag.h"
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <ssht.h>
#include <assert.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))


double eval_laguerre(double z, int n, int alpha)
{
	int k;
	double p1, p2, p3;
	p1 = 1.0 ;
	p2 = 0.0;
	for (k = 1; k <= n; k++){
		p3 = p2;
		p2 = p1;
		p1 = ( (alpha + 2.0 * k - 1.0 - z) * p2 - (alpha + k - 1.0) * p3 ) / k;
	}
	return(p1);
}

double eval_laguerre_rescaled(double z, int n, int alpha, double normfac)
{
	int k;
	double p1, p2, p3;
	p1 = 1.0 / normfac ;
	p2 = 0.0;
	for (k = 1; k <= n; k++){
		p3 = p2;
		p2 = p1;
		p1 = ( (alpha + 2.0 * k - 1.0 - z) * p2 - (alpha + k - 1.0) * p3 ) / k;
	}
	return(p1);
}

long factorial(int n)
{
	int i, res = 1;
	for(i = 1; i <= n; i++)
		res *= i;
	return(res);
}

long factorial_range(int min, int max)
{
	int i, res = 1;
	for(i = min; i <= max; i++)
		res *= i;
	return(res);
}

void flag_spherlaguerre_quadrature(double *roots, double *weights, int N, int alpha)
{
	double infbound ;
	double supbound ;
	const int NITERMAX = 2000;
	const int MAXIT = 250;
	int niter, i, n;
	int facalpha = factorial_range(N+1, N+alpha);

	double h = 1.0 / (double) N;
	double vinf, vsup, p1, p2, pp, z = 0.0, z1, temp;
	double normfac = 1.0;

	// First low-bound
	infbound = h;

	for (n = 0; n < N; n++)
	{
		if( n > 100 )
			h = 0.1;
		supbound = infbound;
		normfac = eval_laguerre_rescaled(z, n, alpha, normfac);

		temp = eval_laguerre_rescaled(infbound, N, alpha, normfac);
		vinf = temp;
		vsup = temp;

		niter = 0;
		while( vinf * vsup >= 0 && niter < NITERMAX ){
			supbound += h;
			vsup = eval_laguerre_rescaled(supbound, N, alpha, normfac);
			niter++;
		}

		niter = 0;
		while( vinf * vsup < 0 && niter < NITERMAX ){
			infbound += h;
			vinf = eval_laguerre_rescaled(infbound, N, alpha, normfac);
			niter++;
		}
		infbound -= h;
		vinf = eval_laguerre_rescaled(infbound, N, alpha, normfac);

		// Linear interpolation for root first estimation
		z = infbound - vinf * (supbound-infbound) / (vsup-vinf);

		// Prepare for next iteration
		infbound = supbound;

		//printf(" %i %f ", n, z);
		//printf(" %1.1e ",normfac);
		// Refine using Newton scheme
		for (i = 1; i <= MAXIT; i++)
		{
			p1 = eval_laguerre_rescaled(z, N, alpha, normfac);
			p2 = eval_laguerre_rescaled(z, N-1, alpha, normfac);
			// Derivative
			pp = (N * p1 - N * p2) / z;
			z1 = z;
			// Newton step
			z = z1 - p1/pp;
			//printf(" %f ", z);
			if( cabs(z - z1) < 1e-16 ){
				//printf("MAXIT = %i err = %2.2e\n",i,cabs(taubis-tau));
				break;
			}
		}
		//printf("\n");
		// Store final root estimate
		roots[n] = z;
		weights[n] = (facalpha / pow(N+1, 2.0)) * z * pow(  (exp(z / 4.0) / normfac) * (exp(z / 4.0) / eval_laguerre_rescaled(z, N+1, alpha, normfac)) , 2.0);
		// Correct weights for classical Gauss-Laguerre quadrature with generalised polynomrials :
		// weights[n] =  (facalpha / pow(N+1, 2.0)) * z / pow( eval_laguerre(z, N+1, alpha), 2.0) ;
	}

	// Specific unstable solution
	if( N == 2 && alpha == 2){
		roots[0] = 2.0;
		roots[1] = 6.0;
		weights[0] = 2.770896037098993835 * pow(roots[0],2.0);
		weights[1] = 5.603177687399098925 * pow(roots[1],2.0);
	}

	// Check order
	for (n = 1; n < N-1; n++)
		if( roots[n] <= roots[n-1] )
			printf("Problem with %ith root! : %f < %f < %f\n",
				n, roots[n-1], roots[n], roots[n+1]);

}



double flag_spherlaguerre_tau(double R, int N)
{
	assert(R > 0.0);
	double Rmax = flag_spherlaguerre_Rmax(N);
	return R / Rmax;
}

double flag_spherlaguerre_Rmax(int N)
{
	assert(N > 1);
	const int alpha = ALPHA;
	double R;

	/*
	double *roots = (double*)calloc(N+1, sizeof(double));
	double *weights = (double*)calloc(N+1, sizeof(double));
	flag_spherlaguerre_quadrature(roots, weights, N+1, alpha);
	double tau_bis = roots[N];
	free(roots);
	free(weights);
	*/

	if( N < 5 ){

		double *roots = (double*)calloc(N, sizeof(double));
		double *weights = (double*)calloc(N, sizeof(double));
		flag_spherlaguerre_quadrature(roots, weights, N, alpha);
		R = roots[N-1];
		free(roots);
		free(weights);

	}else{

		double infbound, supbound =  3.95 * N + 10;
	 	double vinf, vsup, p1, p2, pp, taubis;
	 	double h = (double)N/1000;
	 	const int MAXIT = 250;
	 	int i;

		infbound = supbound;
		double normfac = eval_laguerre(infbound, N, alpha);
		double temp = eval_laguerre_rescaled(infbound, N, alpha, normfac);

		vinf = temp;
		vsup = temp;

		while( vinf * vsup >= 0 && infbound > 0 ){
			infbound -= h;
			vinf = eval_laguerre_rescaled(infbound, N, alpha, normfac);
		}

		while( vinf * vsup < 0 && supbound > 0 ){
			supbound -= h;
			vsup = eval_laguerre_rescaled(supbound, N, alpha, normfac);
		}
		supbound += h;
		vsup = eval_laguerre_rescaled(infbound, N-1, alpha, normfac);

		// Linear interpolation for root first estimation
		R = infbound - vinf * (supbound-infbound) / (vsup-vinf);

		for (i = 1; i <= MAXIT; i++)
		{
			p1 = eval_laguerre_rescaled(R, N, alpha, normfac);
			p2 = eval_laguerre_rescaled(R, N-1, alpha, normfac);
			// Derivative
			pp = (N * p1 - N * p2) / R;
			taubis = R;
			// Newton step
			R = taubis - p1/pp;
			if( cabs(taubis - R) < 1e-16 ){
				//printf("MAXIT = %i err = %2.2e\n",i,cabs(taubis-tau));
				break;
			}
		}
	}

	//printf("R = %4.4e\n",R);
	//printf("Taubis = %4.4e\n",tau_bis);

	return R;
}

void flag_spherlaguerre_sampling(double *nodes, double *weights, double tau, int N)
{
	assert(tau > 0.0);
	assert(N > 1);
	const int alpha = ALPHA;

	flag_spherlaguerre_quadrature(nodes, weights, N, alpha);

	int n;
	for (n=0; n<N; n++){
		//weights[n] *= exp(nodes[n]);// * pow(tau, alpha + 1);
		nodes[n] *= tau;
		//printf("Node %i = %f with weight %f \n",n,nodes[n],weights[n]);
	}

}

void flag_spherlaguerre_allocate_sampling(double **nodes, double **weights, int N)
{
	assert(N > 1);
	*nodes = (double*)calloc(N, sizeof(double));
	*weights = (double*)calloc(N, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
}

void flag_spherlaguerre_analysis(double *fn, const double *f, const double *nodes, const double *weights, double tau, int N)
{
	assert(N > 1);
	int i, n;
	const int alpha = ALPHA;

	double r, factor, lagu0, lagu1, lagu2;

	for(i=0; i<N; i++)
	{
		r = nodes[i]/tau;
		factor = weights[i] * f[i] * exp(-r/2.0) ;

		lagu0 = 0.0;
		lagu1 = 1.0;

		fn[0] += factor * pow(factorial_range(1, alpha), -0.5) * lagu1;

		for (n = 1; n < N; n++)
		{
			lagu2 =
				(
					(alpha + 2 * n - 1 - r) * lagu1 -
					(alpha + n - 1) * lagu0
				) / n;

			fn[n] += factor * pow(factorial_range(n+1, n+alpha), -0.5) * lagu2;

			lagu0 = lagu1;
			lagu1 = lagu2;
		}
	}

}

void flag_spherlaguerre_synthesis(double *f, const double *fn, const double *nodes, int Nnodes, double tau, int N)
{
	flag_spherlaguerre_synthesis_gen(f, fn, nodes, Nnodes, tau, N, ALPHA);
}


void flag_spherlaguerre_synthesis_gen(double *f, const double *fn, const double *nodes, int Nnodes, double tau, int N, int alpha)
{
	assert(N > 1);
	assert(Nnodes > 1);
	int i, n;
	complex double factor, lagu0, lagu1, lagu2, r;

	for (i = 0; i < Nnodes; i++)
	{
		r = nodes[i]/tau;
		factor = exp(-r/2.0);

		lagu0 = 0.0;
		lagu1 = 1.0;

		f[i] += factor * pow(factorial_range(1, alpha), -0.5) * lagu1 * fn[0];

		for (n = 1; n < N; n++)
		{
			lagu2 =
				(
					(alpha + 2 * n - 1 - r) * lagu1 -
					(alpha + n - 1) * lagu0
				) / n;

			f[i] += factor * pow(factorial_range(n+1, n+alpha), -0.5) * lagu2 * fn[n];

			lagu0 = lagu1;
			lagu1 = lagu2;
		}
	}
}


void flag_spherlaguerre_mapped_analysis(complex double *fn, const complex double *f, const double *nodes, const double *weights, double tau, int N, int mapsize)
{
	assert(N > 1);
	assert(mapsize > 1);
	int i, n, l, offset_i, offset_n;
	double factor, lagu0, lagu1, lagu2, r;
	const int alpha = ALPHA;

	double *temp = (double*)calloc(N, sizeof(double));
	//double normfac;

	for(i=0; i<N; i++)
	{

		r = nodes[i]/tau;
		factor = weights[i] * exp(-r/4.0);

		lagu0 = 0.0;
		lagu1 = 1.0 * exp(-r/4.0);

		temp[0] = factor * pow(factorial_range(1, alpha), -0.5) * lagu1;

		for (n = 1; n < N; n++)
		{
			lagu2 =
				(
					(alpha + 2 * n - 1 - r) * lagu1 -
					(alpha + n - 1) * lagu0
				) / n;

			temp[n] = factor * pow(factorial_range(n+1, n+alpha), -0.5) * lagu2;

			lagu0 = lagu1;
			lagu1 = lagu2;
		}

		offset_i = i * mapsize;
		for (n = 0; n < N; n++)
		{
			offset_n = n * mapsize;
			for(l=0; l<mapsize; l++)
			{
				fn[l+offset_n] += f[l+offset_i] * temp[n];
			}
		}
	}

	free(temp);


}

void flag_spherlaguerre_mapped_synthesis(complex double *f, const complex double *fn, const double *nodes, int Nnodes, double tau, int N, int mapsize)
{
	assert(N > 1);
	assert(Nnodes > 1);
	assert(mapsize > 1);
	int i, n, l, offset_n, offset_i;
	const int alpha = ALPHA;

	double r;
	double factor, lagu0, lagu1, lagu2;
	double *temp = (double*)calloc(N, sizeof(double));

	for (i = 0; i < Nnodes; i++)
	{
		r = nodes[i]/tau;
		factor = exp(-r/4.0);
		// was factor = (1.0/r) * exp(-r/2.0) * (1.0/sqrt(tau));

		lagu0 = 0.0;
		lagu1 = 1.0 * exp(-r/4.0);

		temp[0] = factor * pow(factorial_range(1, alpha), -0.5) * lagu1 ;

		for (n = 1; n < N; n++)
		{
			lagu2 =
				(
					(alpha + 2 * n - 1 - r) * lagu1 -
					(alpha + n - 1) * lagu0
				) / n;

			temp[n] = factor * pow(factorial_range(n+1, n+alpha), -0.5) * lagu2;

			lagu0 = lagu1;
			lagu1 = lagu2;

		}

		offset_i = i * mapsize;
		for (n = 0; n < N; n++)
		{
			offset_n = n * mapsize;
			for(l=0; l<mapsize; l++)
			{
				f[l+offset_i] += temp[n] * fn[l+offset_n];
			}
		}
	}

	free(temp);

}


void flag_spherlaguerre_basis(double *KN, const int N, const double *nodes, int Nnodes, double tau)
{
	assert(Nnodes > 0);
	assert(N > 0);
	int i, n;
	const int alpha = ALPHA;
	double factor, lagu0, lagu1, lagu2, r;

	for (i = 0; i < Nnodes; i++)
	{
		r = nodes[i]/tau;
		factor = exp(-r/2.0);

		lagu0 = 0.0;
		lagu1 = 1.0;

		KN[i * N] = factor * pow(factorial_range(1, alpha), -0.5) * lagu1;

		if ( N > 1 ) {

			for (n = 1; n < N; n++)
			{
				lagu2 = ( (alpha + 2 * n - 1 - r) * lagu1 -
						(alpha + n - 1) * lagu0 ) / n;
				lagu0 = lagu1;
				lagu1 = lagu2;
				KN[ i * N + n ] = factor * pow(factorial_range(n+1, n+alpha), -0.5) * lagu2;
			}

		}

	}

}


void flag_spherbessel_sampling(double *nodes, double *weights, double R, int N)
{
	assert(R > 0.0);
	assert(N > 1);
	const int alpha = ALPHA;

	flag_spherlaguerre_quadrature(nodes, weights, N, alpha);
	double tau = R / nodes[N-1];

	int n;
	for (n=0; n<N; n++){
		weights[n] = 1;
		nodes[n] *= tau;
		//printf("Node %i = %f with weight %f \n",n,nodes[n],weights[n]);
	}

}
