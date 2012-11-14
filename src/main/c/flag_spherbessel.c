
// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flag.h"
#include <math.h>
#include <stdlib.h>
#include <complex.h> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>

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
void flag_spherlaguerre2spherbessel(double *flk, const double *fn, double *kvalues, int Nk, int N, int ell, double tau)
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

  free(sbesselslag);

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
void flag_fourierlaguerre2fourierbessel(complex double *flmk, complex double *flmn, double *kvalues, int Nk, int N, int L, double tau)
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
	int k, n, j;
  double weight, hold, result, val, temp, y;

  double *mulk_coefs  = (double*)calloc(N+2, sizeof(double));
  double *vals  = (double*)calloc(N, sizeof(double));

  for(k = 0; k < Nk; k++){
    flag_mulk(mulk_coefs, N+2, ell, kvalues[k], tau);
    for(n = 0; n < N; n++){
      printf("\n");
      for(j = 0; j <= n; j++){
        if( j == 0 ){
          weight =  4 * gsl_sf_fact(ell + 2) * (n + 1.0) * (n + 2.0) / 2.0; // 
        }else{
          weight = - 2 * (ell + 2 + j) * (n - j + 1) * weight / (j * (j + 2)); //  
        }
        //sbesselslag[n * Nk + k] += pow((n+1)*(n+2), -0.5) * weight * mulk_coefs[j+2] ;
        vals[j] = pow((n+1)*(n+2), -0.5) * weight * mulk_coefs[j+2];
      }
      result = 0.0, hold = 0.0, temp = 0.0, y = 0.0;
      for(j = 0; j <= n/2; j++){
        if( 2*j+1 <= n ){
          y = (vals[2*j] + vals[2*j+1]) - hold;
          temp = result + y;
          hold = (temp - result) - y;
          result = temp;
        } else {
          y = vals[2*j] - hold;
          temp = result + y;
          hold = (temp - result) - y;
          result = temp;
        }
        // printf("n = %i, j = %i, val = %3.2e, val = %3.2e, result = ", n, j, vals[2*j], vals[2*j+1]);
        // printf("%3.2e\n",result);
      }
      /*
      bubbleSort(vals, n+1);
      for(j = 0; j <= n; j++)
        printf(" %3.2e \n", vals[j]);
      result = 0.0, hold = 0.0, temp = 0.0, y = 0.0;
      for(j = 0; j <= n/2; j++){
        y = (vals[j] + vals[n-j]) - hold;
        temp = result + y;
        hold = (temp - result) - y;
        result = temp;
        printf("n = %i, j = %i, val = %3.2e, val = %3.2e, result = ", n, j, vals[j], vals[n-j]);
        printf("%3.2e\n",result);
      }
      if( n % 2 != 1 ){
        y = - vals[n/2] - hold;
        temp = result + y;
        hold = (temp - result) - y;
        result = temp;
      }
      */
      sbesselslag[n * Nk + k] = result ;
    }
  }

  free(mulk_coefs);

}

/*!
 * TODO TODO TODO TODO TODO TODO 
 */
void flag_mulk(double *mulk, int n, int ell, double k, double tau)
{
	double PI = 3.141592653589793;
  int j;
  double a, b, fac1, fac2, fac3, rec0, rec1, rec2;
  double c = (double)ell + 1.5;
  double ktilde = tau * k;
  double d = - 4.0 * ktilde * ktilde;
  double fac = pow(tau, 3.0) * sqrt(PI) * pow(ktilde, ell);

  // Even
  a = (ell + 1.0) / 2.0;
  b = ell / 2.0 + 1.0;
  for(j = 0; j <= n/2; j++){
    //printf("j1=%i, 2F1(%f, %f, %f, %f) = ",2*j,a+j,b+j,c,d);fflush(NULL);
    if( j == 0 ){
      rec1 = pow(1-d, -b) * gsl_sf_hyperg_2F1_renorm(b, c-a, c, d/(d-1)) ; //* gsl_sf_fact(2*j+ell) * pow(2, 2*j);
      // rec1 = ( atan(sqrt(-d)) / sqrt(-d) ) / gsl_sf_gamma(c);
      rec0 = rec1;
      //printf("%3.3e\n",rec0);fflush(NULL);
    }else if( j == 1 ){
      rec0 = pow(1-d, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, d/(d-1)) ; //* gsl_sf_fact(2*j+ell) * pow(2, 2*j);
      // rec0 = pow(d-1,-2.0) / gsl_sf_gamma(c);
      //printf("%3.3e\n",rec0);fflush(NULL);
    }else{
      rec2 = rec1;
      rec1 = rec0;
      fac1 = (c-a-j+1)*(c-b-j+1)*(c-a-b-2*j+1) * rec2 ; //* 16*(2*j+ell)*(2*j+ell-1)*(2*j+ell-2)*(2*j+ell-3);
      fac2 = (c-a-b-2*j+2)*( c*(a+b-c+2*j-2) + c - 2.0*(a+j-1)*(b+j-1) 
            + d*( (a+b+2*j-2)*(c-a-b-2*j+2) + 2.0*(a+j-1)*(b+j-1) - c + 1 ) ) * rec1 ; //* 4*(2*j+ell)*(2*j+ell-1);
      fac3 = ( (a+j-1)*(b+j-1)*(c-a-b-2*j+3)*(1.0-d)*(1.0-d) );
      rec0 = -1.0 * (fac1 + fac2) / fac3 ;
      //printf("%3.3e\n",rec0);fflush(NULL);
      //printf("fac1=%5.4f, fac2=%5.4f, fac3=%5.4f\n",fac1,fac2,fac3);fflush(NULL);
    }
    mulk[2*j] = fac * rec0; //* gsl_sf_gamma(2*j+ell+1) * pow(2, 2*j)
  }

  // Odd
  a = ell / 2.0 + 1.0;
  b = (ell + 3) / 2.0;
  for(j = 0; j <= n/2-1; j++){
    //printf("j2=%i, 2F1(%f, %f, %f, %f) = ",2*j+1,a+j,b+j,c,d);fflush(NULL);
    if( j == 0 ){
      rec1 = pow(1-d, -b) * gsl_sf_hyperg_2F1_renorm(b, c-a, c, d/(d-1)) ; //* gsl_sf_fact(2*j+ell+1) * pow(2, 2*j+1);
      //rec1 = pow(1-d, -1.0) / gsl_sf_gamma(c) ;
      rec0 = rec1;
      //printf("%3.3e\n",rec0);fflush(NULL);
    }else if( j == 1 ){
      rec0 = pow(1-d, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, d/(d-1)) ; //* gsl_sf_fact(2*j+ell+1) * pow(2, 2*j+1);
      //rec0 = ((d+3)/(3*pow(1-d, 3.0))) / gsl_sf_gamma(c);//
      //printf("%3.3e\n",rec0);fflush(NULL);
    }else{
      rec2 = rec1;
      rec1 = rec0;
      fac1 = (c-a-j+1)*(c-b-j+1)*(c-a-b-2*j+1) * rec2 ; //* 16*(2*j+ell)*(2*j+ell-1)*(2*j+ell-2)*(2*j+ell-3);
      fac2 = (c-a-b-2*j+2)*( c*(a+b-c+2*j-2) + c - 2.0*(a+j-1)*(b+j-1) 
            + d*( (a+b+2*j-2)*(c-a-b-2*j+2) + 2.0*(a+j-1)*(b+j-1) - c + 1 ) ) * rec1 ; //* 4*(2*j+ell)*(2*j+ell-1);
      fac3 = ( (a+j-1)*(b+j-1)*(c-a-b-2*j+3)*(1.0-d)*(1.0-d) ); 
      rec0 = -1.0 * (fac1 + fac2) / fac3 ;
      //printf("%3.3e\n",rec0);fflush(NULL);
      //printf("fac1=%5.4f, fac2=%5.4f, fac3=%5.4f\n",fac1,fac2,fac3);fflush(NULL);
    }
    mulk[2*j+1] = fac * rec0;//* gsl_sf_gamma(2*j+ell+2) * pow(2, 2*j+1)
  }

  for(j = 0; j <= n; j++){
    //printf("j=%i, mulk(%1.1f, %1.1f, %1.1f, %1.1f) = %6.5e\n",j,(ell+j+1)/2.0,(ell+j)/2.0+1,c,d,mulk[j]);fflush(NULL);
  }

}


void bubbleSort(double numbers[], int array_size)
{
  int i, j;
  double temp;
 
  for (i = (array_size - 1); i > 0; i--)
  {
    for (j = 1; j <= i; j++)
    {
      if ( numbers[j-1] > numbers[j] )
      {
        temp = numbers[j-1];
        numbers[j-1] = numbers[j];
        numbers[j] = temp;
      }
    }
  }
}
