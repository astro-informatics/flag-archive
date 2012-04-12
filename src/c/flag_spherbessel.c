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

} 

void flag_spherbessel_approx(double *flk, const double *f, double *kvalues, int Nk, double *nodes, int Nnodes, int ell)
{
  double deltar;
	int k, r;

	for(k = 0; k < Nk; k++){

		flk[k] = 0.0;
		for(r = 0; r < Nnodes; r++){

      if(r == 0)
        deltar = nodes[1]-nodes[0];
      else if(r == Nnodes-1)
        deltar = nodes[Nnodes-1]-nodes[Nnodes-2];
      else
        deltar = nodes[r+1]-nodes[r];

			flk[k] += pow(r, 2.0) * j_ell( kvalues[k] * nodes[r], ell ) * deltar ;
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

	//result = tau * sqrt(PI) * pow(2, j) * pow(ktilde, ell) 
	//	* gsl_sf_gamma(j + ell + 1) * gsl_sf_gamma(ell + 1.5) 
	//	* gsl_sf_hyperg_2F1(a, b, c, d);

	return result;
}


double j_ell(double X, int l)
{
    int L;
    double JL, AX, AX2;
    double LN2 = 0.6931471805599453094;
    double ONEMLN2 = 0.30685281944005469058277;
    double PID2 = 1.5707963267948966192313217;
    double PID4 = 0.78539816339744830961566084582;
    double ROOTPI12 = 21.269446210866192327578;
    double GAMMA1 = 2.6789385347077476336556;
    double GAMMA2 = 1.3541179394264004169452;
    double PI = 3.141592653589793238463;
    double NU, NU2, BETA, BETA2, COSB;
    double SX, SX2;
    double COTB, COT3B, COT6B, SECB, SEC2B;
    double TRIGARG, EXPTERM, L3;

    AX = abs(X);
    AX2 = pow(AX, 2.0);

    switch(l){
      case 0 :
          if( AX < 0.1 )
             JL = 1.0 - AX2 / 6.0 * (1.0 - AX2 / 20.0);
          else
             JL = sin(AX) / AX;
          break;
      case 1 : 
          if( AX < 0.2 )
             JL = AX / 3.0 * (1.0 - AX2 / 10.0 * ( 1.0 - AX2 / 28.0 ) );
          else
             JL = (sin(AX) / AX - cos(AX)) / AX;
          break;
      case 2 : 
          if( AX < 0.3 )
             JL = AX2 / 15.0 * (1.0 - AX2 / 14.0 * (1.0 - AX2 / 36.0));
          else
             JL=(-3.00 * cos(AX) / AX - sin(AX) * (1.0 - 3.0 / AX2)) / AX;
          break;
      case 3 :
          if( AX < 0.4 )
             JL = AX * AX2 / 105.0 * (1.0 - AX2 / 18.0 * (1.0 - AX2 / 44.0));
          else
             JL = (cos(AX) * (1.0 - 15.0 / AX2) - sin(AX) * (6.0 - 15.0 / AX2) / AX) / AX;
          break;
      case 4 :
          if( AX < 0.6 )
             JL = pow(AX2 ,2.0) / 945.0 * (1.0 - AX2 / 22.0 * (1.0 - AX2 / 52.0));
          else
             JL = (sin(AX) * (1.0 - (45.0 - 105.0 / AX2) / AX2) + cos(AX) * (10.0 - 105.0 / AX2) / AX) / AX;
          break;
      case 5 :
          if( AX < 1.0 )
             JL = pow(AX2, 2.0) * AX / 10395.0 * (1.0 - AX2 / 26.0 * (1.0 - AX2 / 60.0));
          else
             JL = (sin(AX) * (15.0 - (420.0 - 945.0 / AX2) / AX2) / AX - cos(AX) * (1.0 - (105.0 - 945.0 / AX2) / AX2)) / AX;
          break;
      case 6 :
          if( AX < 1.0 )
            JL = pow(AX2, 3.0) / 135135.0 * (1.0 - AX2 / 30.0 * (1.0 - AX2 / 68.0));
          else
             JL = (sin(AX) * ( - 1.0 + (210.0 - (4725.0 - 10395.0 / AX2) / AX2) / AX2)  +  
                  cos(AX) * ( - 21.0 + (1260.0 - 10395.0 / AX2) / AX2) / AX) / AX;
          break;
      default :
    
       NU = 0.50 + L;
       NU2 = pow(NU, 2.0);
       if(AX < 1.0 - 40)
          JL = 0.0;
       else if( (AX2 / L) < 0.5){
          JL = exp(L * log(AX / NU) - LN2 + NU * ONEMLN2 - (1.0 - (1.0 - 3.50 / NU2) / NU2 / 30.0) / 12.0 / NU)  
                / NU * (1.0 - AX2 / (4.0 * NU + 4.0) * (1.0 - AX2 / (8.0 * NU + 16.0) * (1.0 - AX2 / (12.0 * NU + 36.0))));
       }else if((pow(l, 2.0) / AX) < 5.0 - 1){
          BETA = AX - PID2 * (L + 1);
          JL = (cos(BETA) * (1.0 - (NU2 - 0.250) * (NU2 - 2.250) / 8.0 / AX2 * (1.0 - (NU2 - 6.25) * (NU2 - 12.250) / 48.0 / AX2))  
                - sin(BETA) * (NU2 - 0.250) / 2.0 / AX *  (1.0 - (NU2 - 2.250) * (NU2 - 6.250) / 24.0 / AX2 * (1.0 - (NU2 - 12.25) *   
               (NU2 - 20.25) / 80.0 / AX2)) ) / AX   ;
       }else{
          L3 = pow(NU, 0.325);
          if(AX  <  NU - 1.31 * L3){
             COSB = NU / AX;
             SX  =  sqrt(NU2 - AX2);
             COTB = NU / SX;
             SECB = AX / NU;
             BETA = log(COSB + SX / AX);
             COT3B = pow(COTB, 3.0);
             COT6B = pow(COT3B, 2.0);
             SEC2B = pow(SECB, 2.0);
             EXPTERM = ( (2.0 + 3.0 * SEC2B) * COT3B / 24.0  
                   -  ( (4.0 + SEC2B) * SEC2B * COT6B / 16.0  
                   +  ((16.0 - (1512.0 + (3654.0 + 375.0 * SEC2B) * SEC2B) * SEC2B) * COT3B / 5760.0  
                   +  (32.0 + (288.0 + (232.0 + 13.0 * SEC2B) * SEC2B) * SEC2B) * SEC2B * COT6B / 128.0 / NU) * COT6B / NU)  
                   / NU) / NU;
             JL = sqrt(COTB * COSB) / (2.0 * NU) * exp( - NU * BETA + NU / COTB - EXPTERM);
          }else if (AX  >  NU + 1.48 * L3){
             COSB = NU / AX;
             SX = sqrt(AX2 - NU2);
             COTB = NU / SX;
             SECB = AX / NU;
             BETA = acos(COSB);
             COT3B = pow(COTB, 3.0);
             COT6B = pow(COT3B, 2.0);
             SEC2B = pow(SECB, 2.0);
             TRIGARG = NU / COTB - NU * BETA - PID4  
                   - ((2.0 + 3.0 * SEC2B) * COT3B / 24.0   
                   + (16.0 - (1512.0 + (3654.0 + 375.0 * SEC2B) * SEC2B) * SEC2B) * COT3B * COT6B / 5760.0 / NU2) / NU;
             EXPTERM = ( (4.0 + SEC2B) * SEC2B * COT6B / 16.0  
                   - (32.0 + (288.0 + (232.0 + 13.0 * SEC2B) * SEC2B) * SEC2B) * SEC2B *  pow(COT6B, 2.0) / 128.0 / NU2) / NU2;
             JL = sqrt(COTB * COSB) / NU * exp( - EXPTERM) * cos(TRIGARG);
          }else{
             BETA = AX - NU;
             BETA2 = pow(BETA, 2.0);
             SX = 6.0 / AX;
             SX2 = pow(SX, 2.0);
             SECB = pow(SX, 0.3333333333333333);
             SEC2B = pow(SECB, 2.0);
             JL = ( GAMMA1 * SECB  +  BETA * GAMMA2 * SEC2B  
                   - (BETA2 / 18.0 - 1.0 / 45.0) * BETA * SX * SECB * GAMMA1  
                   - ((BETA2 - 1.0) * BETA2 / 36.0 + 1.0 / 420.0) * SX * SEC2B * GAMMA2    
                   + (((BETA2 / 1620.0 - 7.0 / 3240.0) * BETA2 + 1.0 / 648.0) * BETA2 - 1.0 / 8100.0) * SX2 * SECB * GAMMA1  
                   + (((BETA2 / 4536.0 - 1.0 / 810.0) * BETA2 + 19.0 / 11340.0) * BETA2 - 13.0 / 28350.0) * BETA * SX2 * SEC2B * GAMMA2  
                   - ((((BETA2 / 349920.0 - 1.0 / 29160.0) * BETA2 + 71.0 / 583200.0) * BETA2 - 121.0 / 874800.0) *   
                  BETA2 + 7939.0 / 224532000.0) * BETA * SX2 * SX * SECB * GAMMA1) * sqrt(SX) / ROOTPI12;
          }
        }
      break;
    }
        
    if( (X < 0) && (L % 2 != 0) ) 
      JL = -JL;

    return(JL);
}

