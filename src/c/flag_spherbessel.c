// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flag.h"
#include <gsl/gsl_math.h>

/*!
 * Compute the ell-th order spherical Bessel transform from the SLAG transform.
 *
 * \param[out]  flk Output transform in spherical Bessel space.
 * \param[in]  fn Spherical Laguerre transform of some field f.
 * \param[in]  kvalues k scales on which the transform is performed.
 * \param[in]  Nk Number of scales in the array kvalues.
 * \param[in]  N Harmonic band-limit.
 * \param[in]  ell ell order of the spherical Bessel function.
 * \param[in]  tau SLAG rescaling factor.
 * \retval none
 */
void flag_spherlaguerre2spherbessel(double *flk, double *fn, double *kvalues, int Nk, int N, int ell, double tau)
{
	double PI = 3.141592653589793;
	double *sbesselslag  = (double*)calloc(N*Nk, sizeof(double));

	flag_sbesselslag(sbesselslag, ell, kvalues, Nk, N, tau);

	int k, n;
	for(k = 0; k < Nk; k++){
		for(n = 0; n < N; n++){
			flk[k] += sqrt(2.0/PI) * fn[n] * sbesselslag[n * Nk + k];
		}
	}

} 

/*!
 * Compute the Fourier-Bessel transform from the Fourier-Laguerre transform.
 *
 * \param[out]  flmk Output Fourier-Bessel transform.
 * \param[in]  flmn Spherical Laguerre transform of some 3D field f.
 * \param[in]  kvalues k scales on which the transform is performed.
 * \param[in]  Nk Number of scales in the array kvalues.
 * \param[in]  N Radial harmonic band-limit.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  tau SLAG rescaling factor.
 * \retval none
 */
void flag_FLAG2FBESS(complex double *flmk, complex double *flmn, double *kvalues, int Nk, int N, int L, double tau)
{

} 

/*!
 * Compute SLAG transform of spherical Bessel function j_ell.
 *
 * \param[out]  sbesselslag Output transform in SLAG space.
 * \param[in]  ell ell order of the spherical Bessel function.
 * \param[in]  kvalues k scales on which the transform is performed.
 * \param[in]  Nk Number of scales in the array kvalues.
 * \param[in]  N Harmonic band-limit.
 * \param[in]  tau SLAG rescaling factor.
 * \retval none
 */
void flag_sbesselslag(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau)
{
	double *weights = (double*)calloc( N*(N+1)/2 , sizeof(double));

	int k, n, j;
	for(k = 0; k < Nk; k++){
		for(n = 0; n < N; n++){
			printf("\n n = %i \n", n);
			for(j = n; j > 0; j--)
				weights[j] = - n * weights[j-1] / (j*j);
			weights[0] = 1;
			for(j = 0; j <= n; j++){
				printf(" %f ",weights[j]);
				sbesselslag[n * Nk + k] += weights[n] * flag_mujlk(j + 1, ell, kvalues[k], tau) ;
			}
		}
	}
	free(weights);

}

/*!
 * Compute moments of the weighted spherical Bessel functions.
 *
 * \param[in]  j j-th moment.
 * \param[in]  ell ell order of the spherical Bessel function.
 * \param[in]  k scale on which the transform is performed.
 * \param[in]  tau SLAG rescaling factor.
 * \retval flag_mujlk the final number.
 */
double flag_mujlk(int j, int ell, double k, double tau)
{
	double result;
	double PI = 3.141592653589793;

	double a = (j + ell + 1) / 2.0;
	double b = 1.0 + (j + ell) / 2.0;
	double c = (double)ell + 1.5;
	double ktilde = tau * k;
	double d = -4.0 * ktilde * ktilde;

	result = tau * sqrt(PI) * pow(2, j) * pow(ktilde, ell) 
		* gsl_sf_gamma(j + ell + 1) * gsl_sf_gamma(ell + 1.5) 
		* gsl_sf_hyperg_2F1(a, b, c, d);

	return result;
}