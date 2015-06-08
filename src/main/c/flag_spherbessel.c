
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
  int k, n, j, m;
  long double hold, result, temp, y, fac, corrfac, f_jl2bis;
  long double cw_jp, c_jp, w_jl, ktilde, a, b, c, z, f_jl2_even, f_jl2_odd, f_jl1_even, f_jl1_odd, f_jl0, f_jl1, f_jl2, T0, T1, T2;
  long double PI = 3.141592653589793;

  double *mulk_coefs  = (double*)calloc(N+2, sizeof(double));
  long double *vals  = (double*)calloc(N, sizeof(long double));

  for(k = 0; k < Nk; k++){

    //flag_mulk(mulk_coefs, N+2, ell, kvalues[k], tau);
  
    ktilde = tau * kvalues[k];

    for(n = 0; n < N; n++){


      for(j = 0; j <= n; j++){

	if ( (j % 2) == 0 ){ // j even
	  m = (long double)j / 2.0;
	  a = ((long double)ell + 1) / 2.0;
	  b = (long double)ell / 2.0 + 1;
	  c = (long double)ell + 1.5;
	  z = - 4.0 * pow(ktilde, 2.0);
	  f_jl0 = f_jl1_even;
	  f_jl1 = f_jl2_even;
	  //printf("\n j=%i even ; a=%f, b=%f, c=%f",j,a,b,c);
	} else { // j odd
	  m = ((long double)j - 1) / 2.0;
	  a = (long double)ell / 2.0 + 1.0;
	  b = ((long double)ell + 3) / 2.0;
	  c = (long double)ell + 1.5;
	  z = - (long double)4 * pow(ktilde, 2.0);
	  f_jl0 = f_jl1_odd;
	  f_jl1 = f_jl2_odd;
	  //printf("\n j=%i odd ; a=%f, b=%f, c=%f",j,a,b,c);
	}

        if( j == 0 ){
	  fac = 1.0;
          c_jp = ((long double)n + 2) * (n + 1) / 2.0 ;
          w_jl = sqrt(PI) * pow(ktilde, ell) * pow(tau, 3.0) * gsl_sf_fact(ell + 2) ; 
	  cw_jp = c_jp * w_jl;
	  corrfac = 1.0;//pow(1-z, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, z/(z-1));
	  f_jl1_even = pow(1-z,-a) * gsl_sf_hyperg_2F1_renorm(c-b,a,c,z/(z-1)) / corrfac; 
	  f_jl2_even = pow(1-z, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, z/(z-1)) / corrfac; // pow(1-z,-a-1) * gsl_sf_hyperg_2F1_renorm(c-b-1,a+1,c,z/(z-1));
	  f_jl2 = f_jl2_even;
	  sbesselslag[n * Nk + k] += pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+2) * fac * cw_jp * f_jl2 ;
	  printf("k=%f n=%i j=%i difference=%Le\n", kvalues[k], n, j,  pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+2) * fac * cw_jp * f_jl2  );

        }else{
	  fac = ((long double)ell + 2 + j) * (-1) * ((long double)n - j + 1.0) / ((long double)j * (j + 2.0));
	  if( j == 1 ){
	    f_jl2_odd = pow(1-z, -b) * gsl_sf_hyperg_2F1_renorm(b, c-a, c, z/(z-1)) / corrfac; //
	    f_jl1_odd = pow(1-z, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, z/(z-1)) / corrfac; 
	    f_jl2bis = f_jl1_odd;
	    sbesselslag[n * Nk + k] += pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+1) * cw_jp * (f_jl2 + 2 * fac * f_jl2bis);
	    printf("k=%f n=%i j=%i difference=%Le\n", kvalues[k], n, j,  pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+1) * cw_jp * (f_jl2 + 2 * fac * f_jl2bis) );
	    f_jl2 = f_jl2bis;
	  } else {
	    T0 = (c - a - m) * (c - b - m) * (c - a - b - 2*m - 1.0);
	    T1 = (c-a-b-2*m)*( c*(a+b-c+2*m) + c - 2*(a+m)*(b+m) + z*( (a+b+2*m)*(c-a-b-2*m) + 2*(a+m)*(b+m) - c + 1.0 ) );
	    T2 = (a + m) * (b + m) * (c - a - b - 2*m + 1.0) * pow(1 - z, 2.0) ;
	    if ( (j % 2) == 0 ){ // j even
	      f_jl2 = (-T1/T2) * f_jl1 + (-T0/T2) * f_jl0 ;
	      if( j == n ){
		sbesselslag[n * Nk + k] += pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+2) * fac * cw_jp * f_jl2 ;
		printf("k=%f n=%i j=%i difference=%Le\n", kvalues[k], n, j,  pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+2) * fac * cw_jp * f_jl2  );
	      }
	    }else{
	      f_jl2bis = ( (-T1/T2) * f_jl1 + (-T0/T2) * f_jl0 ) ;
	      sbesselslag[n * Nk + k] += pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+1) * cw_jp * (f_jl2 + 2 * fac * f_jl2bis);
	      printf("k=%f n=%i j=%i difference=%Le\n", kvalues[k], n, j,  pow((n+1)*(n+2), -0.5) * pow(2,j+1) * cw_jp *  (f_jl2 + 2 * fac * f_jl2bis) );
	      f_jl2 = f_jl2bis;
	    }
	  }
	  cw_jp = fac * cw_jp ;
          c_jp = - ((long double)n - j + 1.0) / ((long double)j * (j + 2.0)) * c_jp ; 
	  w_jl = ((long double)ell + 2 + j) * w_jl;
	}

	if ( (j % 2) == 0 ){ // j even
	  f_jl1_even = f_jl2_even;
	  f_jl2_even = f_jl2;
	} else { // j odd
	  f_jl1_odd = f_jl2_odd;
	  f_jl2_odd = f_jl2;
	}

	//vals[j] = pow((n+1)*(n+2), -0.5) * cw_jp * ( pow(2,j+2) * f_jl2 );
	//printf("\n  cw_jp=%e  f_jl2=%e  val=%e",  cw_jp, f_jl2, vals[j]);
	//printf("\n  k=%f  n=%i  j=%i  cw_jp=%Le  f_jl2=%Le  val=%Le", kvalues[k], n, j, cw_jp, f_jl2, vals[j]);
	//sbesselslag[n * Nk + k] += vals[j]  ;

      }    
      printf("sbesselslag[n * Nk + k] = %e\n", sbesselslag[n * Nk + k]);
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
    mulk[2*j] =  rec0; //fac * //* gsl_sf_gamma(2*j+ell+1) * pow(2, 2*j)
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
    mulk[2*j+1] = rec0;//fac * * gsl_sf_gamma(2*j+ell+2) * pow(2, 2*j+1)
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


void flag_sbesselslag_backup(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau)
{
	int k, n, j;
  double weight, hold, result, temp, y;

  double *mulk_coefs  = (double*)calloc(N+2, sizeof(double));
  double *vals  = (double*)calloc(N, sizeof(double));

  for(k = 0; k < Nk; k++){
    flag_mulk(mulk_coefs, N+2, ell, kvalues[k], tau);
    for(n = 0; n < N; n++){
      printf("\n");
      for(j = 0; j <= n; j++){
        if( j == 0 ){
          weight =  pow(2, 2) * gsl_sf_fact(ell + 2) * (n + 1.0) * (n + 2.0) / 2.0 ; // 
        }else{
          weight = - 2 * (ell + 2 + j) * (n - j + 1) * weight / (j * (j + 2)); //  
        }
        sbesselslag[n * Nk + k] += pow((n+1)*(n+2), -0.5) * weight * mulk_coefs[j+2] ;
        printf(" %3.0e ",weight);
        vals[j] = pow((n+1)*(n+2), -0.5) * ( weight * mulk_coefs[j+2] );
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
void flag_sbesselslag_backup2(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau)
{
  int k, n, j, m;
  double hold, result, temp, y;
  double c_jp, w_jl, ktilde, a, b, c, z, f_jl2_even, f_jl2_odd, f_jl1_even, f_jl1_odd, f_jl0, f_jl1, f_jl2, T0, T1, T2;
  double PI = 3.141592653589793;

  double *mulk_coefs  = (double*)calloc(N+2, sizeof(double));
  double *vals  = (double*)calloc(N, sizeof(double));

  for(k = 0; k < Nk; k++){

    flag_mulk(mulk_coefs, N+2, ell, kvalues[k], tau);
  
    ktilde = tau * kvalues[k];

    for(n = 0; n < N; n++){

      printf("\n n = %i", n);

      for(j = 0; j <= n; j++){

	if ( (j % 2) == 0 ){ // j even
	  m = j / 2.0;
	  a = (ell + 1) / 2.0;
	  b = ell / 2.0 + 1;
	  c = ell + 1.5;
	  z = - 4.0 * pow(ktilde, 2.0);
	  f_jl0 = f_jl1_even;
	  f_jl1 = f_jl2_even;
	  //printf("\n j=%i even ; a=%f, b=%f, c=%f",j,a,b,c);
	} else { // j odd
	  m = ((double)j - 1) / 2.0;
	  a = (double)ell / 2.0 + 1.0;
	  b = ((double)ell + 3) / 2.0;
	  c = (double)ell + 1.5;
	  z = - (double)4 * pow(ktilde, 2.0);
	  f_jl0 = f_jl1_odd;
	  f_jl1 = f_jl2_odd;
	  //printf("\n j=%i odd ; a=%f, b=%f, c=%f",j,a,b,c);
	}

        if( j == 0 ){
          c_jp = (n + 2) * (n + 1) / 2.0 ;
          w_jl =  4 * sqrt(PI) * pow(ktilde, ell) * pow(tau, 1.5) * gsl_sf_fact(ell + 2) ; 
	  f_jl1_even = pow(1-z, -b) * gsl_sf_hyperg_2F1_renorm(b, c-a, c, z/(z-1)); // pow(1-z,-a) * gsl_sf_hyperg_2F1_renorm(c-b,a,c,z/(z-1));
	  f_jl2_even = pow(1-z, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, z/(z-1)); // pow(1-z,-a-1) * gsl_sf_hyperg_2F1_renorm(c-b-1,a+1,c,z/(z-1));
	  f_jl2 = f_jl2_even;
        }else{
          c_jp = - (n - j + 1) * c_jp / (j * (j + 2)); 
	  w_jl = 2 * (ell + 2 + j) * w_jl;
	  if( j == 1 ){
	    f_jl2_odd = pow(1-z, -b) * gsl_sf_hyperg_2F1_renorm(b, c-a, c, z/(z-1)); // pow(1-z,-a) * gsl_sf_hyperg_2F1_renorm(c-b,a,c,z/(z-1));
	    f_jl1_odd = pow(1-z, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, z/(z-1)) ; //* gsl_sf_fact(2*j+ell+1) * pow(2, 2*j+1);
	    //f_jl2_odd = pow(1-z, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, z/(z-1)); // pow(1-z,-a-1) * gsl_sf_hyperg_2F1_renorm(c-b-1,a+1,c,z/(z-1));
	    f_jl2 = f_jl1_odd;
	  } else {
	    T0 = (c - a - m) * (c - b - m) * (c - a - b - 2*m - 1.0);
	    T1 = (c-a-b-2*m)*( c*(a+b-c+2*m) + c - 2*(a+m)*(b+m) + z*( (a+b+2*m)*(c-a-b-2*m) + 2*(a+m)*(b+m) - c + 1.0 ) );
	    T2 = (a + m) * (b + m) * (c - a - b - 2*m + 1.0) * pow(1 - z, 2.0) ;
	    f_jl2 = (-T1/T2) * f_jl1 + (-T0/T2) * f_jl0 ;
	  }

	  if ( (j % 2) == 0 ){ // j even
	    f_jl1_even = f_jl2_even;
	    f_jl2_even = f_jl2;
	  } else { // j odd
	    f_jl1_odd = f_jl2_odd;
	    f_jl2_odd = f_jl2;
	  }

	  vals[j] = c_jp * w_jl * f_jl2 ;
	  printf("\n c_jp = %e  w_jl = %e  f_jl2 = %e  mu = %e  val = %f", c_jp, w_jl, f_jl2);
	  sbesselslag[n * Nk + k] += vals[j];
	  sbesselslag[n * Nk + k] *= pow((n+1)*(n+2), -0.5);
        }

      }
      /*
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
      */
      //sbesselslag[n * Nk + k] = result ;
    }
  }

  free(mulk_coefs);

}



void flag_sbesselslag_backup3(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau)
{
  int k, n, j, m;
  long double hold, result, temp, y, fac, corrfac, f_jl2bis;
  long double cw_jp, c_jp, w_jl, ktilde, a, b, c, z, f_jl2_even, f_jl2_odd, f_jl1_even, f_jl1_odd, f_jl0, f_jl1, f_jl2, T0, T1, T2;
  long double PI = 3.141592653589793;

  double *mulk_coefs  = (double*)calloc(N+2, sizeof(double));
  long double *vals  = (double*)calloc(N, sizeof(long double));

  for(k = 0; k < Nk; k++){

    //flag_mulk(mulk_coefs, N+2, ell, kvalues[k], tau);
  
    ktilde = tau * kvalues[k];

    for(n = 0; n < N; n++){

      //printf("\n n = %i", n);

      for(j = 0; j <= n; j++){

	if ( (j % 2) == 0 ){ // j even
	  m = (long double)j / 2.0;
	  a = ((long double)ell + 1) / 2.0;
	  b = (long double)ell / 2.0 + 1;
	  c = (long double)ell + 1.5;
	  z = - 4.0 * pow(ktilde, 2.0);
	  f_jl0 = f_jl1_even;
	  f_jl1 = f_jl2_even;
	  //printf("\n j=%i even ; a=%f, b=%f, c=%f",j,a,b,c);
	} else { // j odd
	  m = ((long double)j - 1) / 2.0;
	  a = (long double)ell / 2.0 + 1.0;
	  b = ((long double)ell + 3) / 2.0;
	  c = (long double)ell + 1.5;
	  z = - (long double)4 * pow(ktilde, 2.0);
	  f_jl0 = f_jl1_odd;
	  f_jl1 = f_jl2_odd;
	  //printf("\n j=%i odd ; a=%f, b=%f, c=%f",j,a,b,c);
	}

        if( j == 0 ){
	  fac = 1.0;
          c_jp = ((long double)n + 2) * (n + 1) / 2.0 ;
          w_jl = sqrt(PI) * pow(ktilde, ell) * pow(tau, 3.0) * gsl_sf_fact(ell + 2) ; 
	  cw_jp = c_jp * w_jl;
	  corrfac = 1.0;//pow(1-z, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, z/(z-1));
	  f_jl1_even = pow(1-z,-a) * gsl_sf_hyperg_2F1_renorm(c-b,a,c,z/(z-1)) / corrfac; 
	  f_jl2_even = pow(1-z, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, z/(z-1)) / corrfac; // pow(1-z,-a-1) * gsl_sf_hyperg_2F1_renorm(c-b-1,a+1,c,z/(z-1));
	  f_jl2 = f_jl2_even;
        }else{
	  fac = ((long double)ell + 2 + j) * (-1) * ((long double)n - j + 1.0) / ((long double)j * (j + 2.0));
	  if( j == 1 ){
	    f_jl2_odd = pow(1-z, -b) * gsl_sf_hyperg_2F1_renorm(b, c-a, c, z/(z-1)) / corrfac; //
	    f_jl1_odd = pow(1-z, -a-1) * gsl_sf_hyperg_2F1_renorm(a+1, c-b-1, c, z/(z-1)) / corrfac; 
	    f_jl2bis = f_jl1_odd;
	    sbesselslag[n * Nk + k] += pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+1) * cw_jp * (f_jl2 + 2 * fac * f_jl2bis);
	    f_jl2 = f_jl2bis;
	  } else {
	    T0 = (c - a - m) * (c - b - m) * (c - a - b - 2*m - 1.0);
	    T1 = (c-a-b-2*m)*( c*(a+b-c+2*m) + c - 2*(a+m)*(b+m) + z*( (a+b+2*m)*(c-a-b-2*m) + 2*(a+m)*(b+m) - c + 1.0 ) );
	    T2 = (a + m) * (b + m) * (c - a - b - 2*m + 1.0) * pow(1 - z, 2.0) ;
	    if ( (j % 2) == 0 ){ // j even
	      f_jl2 = (-T1/T2) * f_jl1 + (-T0/T2) * f_jl0 ;
	      if( j == n ){
		sbesselslag[n * Nk + k] += pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+2) * fac * cw_jp * f_jl2 ;
	      }
	    }else{
	      f_jl2bis = ( (-T1/T2) * f_jl1 + (-T0/T2) * f_jl0 ) ;
	      sbesselslag[n * Nk + k] += pow((n+1)*(n+2), -0.5) * corrfac * pow(2,j+1) * cw_jp * (f_jl2 + 2 * fac * f_jl2bis);
	      printf("n=%i j=%i difference=%Le\n", n, j,  pow((n+1)*(n+2), -0.5) * pow(2,j+1) * cw_jp *  (f_jl2 + 2 * fac * f_jl2bis) );
	      f_jl2 = f_jl2bis;
	    }
	  }
	  cw_jp = fac * cw_jp ;
          c_jp = - ((long double)n - j + 1.0) / ((long double)j * (j + 2.0)) * c_jp ; 
	  w_jl = ((long double)ell + 2 + j) * w_jl;
	}

	if ( (j % 2) == 0 ){ // j even
	  f_jl1_even = f_jl2_even;
	  f_jl2_even = f_jl2;
	} else { // j odd
	  f_jl1_odd = f_jl2_odd;
	  f_jl2_odd = f_jl2;
	}

	//vals[j] = pow((n+1)*(n+2), -0.5) * cw_jp * ( pow(2,j+2) * f_jl2 );
	//printf("\n  cw_jp=%e  f_jl2=%e  val=%e",  cw_jp, f_jl2, vals[j]);
	//printf("\n  k=%f  n=%i  j=%i  cw_jp=%Le  f_jl2=%Le  val=%Le", kvalues[k], n, j, cw_jp, f_jl2, vals[j]);
	//sbesselslag[n * Nk + k] += vals[j]  ;
	

      }    
      printf("sbesselslag[n * Nk + k] = %e\n", sbesselslag[n * Nk + k]);
    }
  }

  free(mulk_coefs);

}
