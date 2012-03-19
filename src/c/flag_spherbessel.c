// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flag.h"
#include <gsl/gsl_math.h>

void flag_spherlaguerre2spherbessel(complex double *flk, complex double *fn, double *kvalues, int Nk, int N, int ell, double tau)
{
	double PI = 3.141592653589793;
	double *sbesselslag;
	int k, n;
	for(k = 0; k < Nk; k++){
		for(n = 0; n < N; n++){
			flk[k] += sqrt(2.0/PI) * fn[n] * sbesselslag[n * Nk + k];
		}
	}

} 

void flag_FLAG2FBESS(complex double *flmk, complex double *flmn, double *k, int Nk, int N, int L, double tau)
{

} 

void flag_sbesselslag(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau)
{
	double *coefs;
	int k, n, j;
	for(k = 0; k < Nk; k++){
		for(n = 0; n < N; n++){
			for(j = 0; j <= n; j++){
				sbesselslag[n * Nk + k] += coefs[n] * flag_mujlk(j + 1, ell, kvalues[k], tau) ;
			}
		}
	}
}

double flag_mujlk(int j, int ell, double k, double tau)
{
	double result;

	double a = (j + ell + 1) / 2.0;
	double b = 1.0 + (j + ell) / 2.0;
	double c = (double)ell + 1.5;
	double ktilde = tau * k;
	double d = -4.0 * ktilde * ktilde;

	//result = tau * ... * gsl_sf_hyperg_2F1(a, b, c, d);

	return result;
}