
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h> 

#include "flag_error.h"

void flag_spherlaguerre_quadrature(double *roots, double *weights, int N){

	int i, n, k;

	const int MAXIT = 10;         
                                          
	//double EPS = 0.00000001;         
	const double C1 = 0.9084064; 
	const double C2 = 0.05214976;
	const double C3 = 0.002579930;
	const double C4 = 0.003986126;
	const double PI = 3.141592653589793;

	const double anu = 4.0*N + 2.0;

	int rhs = 4*N+3;
	double rhsb, r3, r2, theta, z, p1, p2, p3, pp, z1, ppp;

	for (n = 0; n < N; n++)
	{
		rhs = rhs - 4;
		rhsb = rhs * PI/anu;
		r3 = pow(rhsb, 1.0 / 3.0);
		r2 = pow(r3, 2.0);
		theta = r3 * (C1 + r2 * 
			(C2 + r2 * (C3 + r2 * C4)));
		z = anu * pow(cos(theta), 2.0);

		for (i = 1; i <= MAXIT; i++)
		{
			p1 = 1.0;
			p2 = 0.0;
			for (k = 1; k <= N; k++)
			{
			    p3 = p2;
			    p2 = p1;
			    p1 = ((2.0*k-1.0-z)*p2-(k-1.0)*p3)/k;
			}
			pp = (N*p1-N*p2)/z;
			z1 = z;
			z = z1-p1/pp;
			   	/*if( abs(z[n]-z1[n]) < EPS*abs(z[n]) ){
			   		unfinished[n] = 1;
			   		printf("%i %f \n",n,abs(z[n]-z1[n])/z[n]);
			    }*/
			p3 = p2;
			p2 = p1;
			ppp = ((2.0*N+1.0-z)*p2-N*p3)/(N+1.0);
	   } 
		roots[n] = z;
		weights[n] =  exp(z)*z / (pow(N+1,2.0)*pow(ppp,2.0)); 
	}

}

double flag_spherlaguerre_tau(double R, int N){
	
	int k, i;
	int MAXIT = 5;
	const double C1 = 0.9084064; 
	const double C2 = 0.05214976;
	const double C3 = 0.002579930;
	const double C4 = 0.003986126;
	const double PI = 3.141592653589793;

	const double anu = (4.0*N + 2.0);
	const double rhs = 3*PI/anu;

	const double r3 = pow(rhs, 1.0/3.0);
	const double r2 = pow(r3, 2.0);
	const double theta = r3*(C1+r2*(C2+r2*(C3+r2*C4)));
	double tau = anu*pow(cos(theta), 2.0);
	
	double p1, p2, p3, pp, z1;

	for (i=1; i<=MAXIT; i++){
		p1 = 1.0;
		p2 = 0.0;
		for (k=1; k<=N; k++){
			p3 = p2;
			p2 = p1;
			p1 = ((2.0*k-1.0-tau)*p2-(k-1.0)*p3)/k;
		}
		pp = (N*p1-N*p2)/tau;
		z1 = tau;
		tau = z1-p1/pp; 
	}

	return R / tau;
}

void flag_spherlaguerre_sampling(double *nodes, double *weights, double R, int N){
	
	flag_spherlaguerre_quadrature(nodes, weights, N);

	double tau = R / nodes[N-1];
	
	int n;
	for (n=0; n<N; n++){
		nodes[n] *= tau;
		weights[n] *= tau;
	}

}

void flag_spherlaguerre_allocate_sampling(double *nodes, double *weights, int N){
	
	nodes = (double*)calloc(N, sizeof(double));
	weights = (double*)calloc(N, sizeof(double));

}

void flag_spherlaguerre_analysis(double *fn, const double *f, const double *nodes, const double *weights, int N){

	int i, n;

	const double R = nodes[N-1];
	const double tau = flag_spherlaguerre_tau(R, N);

	for(i=0; i<N; i++)
	{

		double r = nodes[i]/tau;
		double factor = f[i] * weights[i] * exp(-r/2.0) * (1.0/sqrt(tau));

		double lagu0 = 1.0;
		double lagu1 = 1.0 - r;
		double lagu2;

		fn[0] += factor * lagu0;
		fn[1] += factor * lagu1;

		for (n = 2; n < N; n++) 
		{ 
			lagu2 = 
				( 
					((2 * n - 1) - r) * lagu1 - 
					(n - 1) * lagu0
				) / n;

			fn[n] += factor * lagu2;

			lagu0 = lagu1;
			lagu1 = lagu2;
		}
	}

	
}

void flag_spherlaguerre_synthesis(double *f, const double *fn, const double *nodes, int N){
	
	int i, n;

	const double R = nodes[N-1];
	const double tau = flag_spherlaguerre_tau(R, N);

	for (i = 0; i < N; i++)
	{

		double r = nodes[i]/tau;
		double factor = exp(-r/2.0) * (1.0/sqrt(tau));

		double lagu0 = 1.0;
		double lagu1 = 1.0 - r;
		double lagu2;

		f[i] += factor * lagu0 * fn[0];
		f[i] += factor * lagu1 * fn[1];

		for (n = 2; n < N; n++) 
		{ 
			lagu2 = 
				( 
					((2 * n - 1) - r) * lagu1 - 
					(n - 1) * lagu0
				) / n;

			f[i] += factor * lagu2 * fn[n];

			lagu0 = lagu1;
			lagu1 = lagu2;
		}
		
	}

}
