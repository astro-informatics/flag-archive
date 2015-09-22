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

int ssht_fr_size_mw(int L)
{// In case we want to extend to various sampling schemes.
	assert(L > 0);
	int mapsize = L*(2*L-1); // MW sampling scheme
	return mapsize;
}

int ssht_flm_size(int L)
{// In case we want to extend to various sampling schemes
	assert(L > 0);
	int mapsize = L*L;
	return mapsize;
}

int flag_core_flmn_size(int L, int N)
{// In case we want to extend to various sampling schemes
	assert(L > 0);
	assert(N > 1);
	return ssht_flm_size(L)*N;
}

int flag_core_f_size_mw(int L, int N)
{// In case we want to extend to various sampling schemes
	assert(L > 0);
	assert(N > 1);
	return ssht_fr_size_mw(L) * (N);
}

void flag_core_allocate_flmn(complex double **flmn, int L, int N)
{
	assert(L > 0);
	assert(N > 1);
	int flmsize = ssht_flm_size(L);
	long totalsize = N*flmsize;
	*flmn = (complex double*)calloc(totalsize, sizeof(complex double));
	assert(flmn != NULL);
}

void flag_core_allocate_f(complex double **f, int L, int N)
{
	assert(L > 0);
	assert(N > 1);
	int frsize = ssht_fr_size_mw(L);
	long totalsize = (N) * frsize;
	*f = (complex double*)calloc(totalsize, sizeof(complex double));
	assert(f != NULL);
}

void flag_core_allocate_f_real(double **f, int L, int N)
{
	assert(L > 0);
	assert(N > 1);
	int frsize = ssht_fr_size_mw(L);
	long totalsize = (N) * frsize;
	*f = (double*)calloc(totalsize, sizeof(double));
	assert(f != NULL);
}

void flag_core_analysis(complex double *flmn,
		const complex double *f,
		int L, double tau, int N, int spin)
{
	assert(L > 0);
	assert(N > 1);
	//const int alpha = ALPHA;
	int verbosity = 0;
	int n;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
	int offset_lm, offset_r;

	complex double *flmr;
	flag_core_allocate_flmn(&flmr, L, N);

	for (n = 0; n < N; n++){
		//printf("> Analysis: layer %i on %i\n", n+1,N+1);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_forward_sov_conv_sym(flmr + offset_lm, f + offset_r, L, spin, dl_method, verbosity);
	}

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, tau, N);
	//printf("> Mapped spherical Laguerre transform...");
	fflush(NULL);
	flag_spherlaguerre_mapped_analysis(flmn, flmr, nodes, weights, tau, N, flmsize);
	//printf("done\n");
	free(nodes);
	free(weights);
    free(flmr);

}

void flag_core_synthesis(complex double *f,
		const complex double *flmn,
		const double *nodes, int Nnodes,
		int L, double tau, int N, int spin)
{
	assert(L > 0);
	assert(N > 1);
	assert(nodes != NULL);
	//const int alpha = ALPHA;
	int verbosity = 0;
	int n, offset_lm, offset_r;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;

	complex double *flmr;
	flag_core_allocate_flmn(&flmr, L, Nnodes);
	//printf("> Mapped spherical Laguerre transform...");fflush(NULL);
	flag_spherlaguerre_mapped_synthesis(flmr, flmn, nodes, Nnodes, tau, N, flmsize);
	//printf("done\n");

	for (n = 0; n < Nnodes; n++){
		//printf("> Synthesis: layer %i on %i\n",n+1,Nnodes);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_inverse_sov_sym(f + offset_r, flmr + offset_lm, L, spin, dl_method, verbosity);
	}

    free(flmr);
}

void flag_core_analysis_real(complex double *flmn,
		const double *f, int L, double tau, int N)
{
	assert(L > 0);
	assert(N > 1);
	int verbosity = 0;
	int n, offset_lm, offset_r;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;

	complex double *flmr;
	flag_core_allocate_flmn(&flmr, L, N);

	for (n = 0; n < N; n++){
		//printf("> Analysis: layer %i on %i\n",n+1,N);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_forward_sov_conv_sym_real(flmr + offset_lm, f + offset_r, L, dl_method, verbosity);
	}

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, tau, N);
	//printf("> Mapped spherical Laguerre transform...");
	fflush(NULL);
	flag_spherlaguerre_mapped_analysis(flmn, flmr, nodes, weights, tau, N, flmsize);
	//printf("done\n");
	free(nodes);
	free(weights);
    free(flmr);
}

void flag_core_synthesis_real(double *f,
		const complex double *flmn,
		const double *nodes, int Nnodes,
		int L, double tau, int N)
{
	assert(L > 0);
	assert(N > 1);
	int verbosity = 0;
	int n, offset_lm, offset_r;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;

	complex double *flmr;
	//printf("> Mapped spherical Laguerre transform...");fflush(NULL);
	flag_core_allocate_flmn(&flmr, L, Nnodes);
	flag_spherlaguerre_mapped_synthesis(flmr, flmn, nodes, Nnodes, tau, N, flmsize);
	//printf("done\n");

	for (n = 0; n < Nnodes; n++){
		//printf("> Synthesis: layer %i on %i\n",n+1,N);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_inverse_sov_sym_real(f + offset_r, flmr + offset_lm, L, dl_method, verbosity);
	}

    free(flmr);
}


double j_ell(double X, int l)
{
    int L = l;
    double JL, AX, AX2;
    double LN2 = 0.6931471805599453094;
    double ONEMLN2 = 0.30685281944005469058277;
    double PID2 = 1.5707963267948966192313217;
    double PID4 = 0.78539816339744830961566084582;
    double ROOTPI12 = 21.269446210866192327578;
    double GAMMA1 = 2.6789385347077476336556;
    double GAMMA2 = 1.3541179394264004169452;
    //double PI = 3.141592653589793238463;
    double NU, NU2, BETA, BETA2, COSB;
    double SX, SX2;
    double COTB, COT3B, COT6B, SECB, SEC2B;
    double TRIGARG, EXPTERM, L3;

    AX = pow(pow(X, 2.0), 0.5);
    AX2 = pow(X, 2.0);

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
       else if( (AX2 / (double)L) < 0.5){
          JL = exp(L * log(AX / NU) - LN2 + NU * ONEMLN2 - (1.0 - (1.0 - 3.50 / NU2) / NU2 / 30.0) / 12.0 / NU)
                / NU * (1.0 - AX2 / (4.0 * NU + 4.0) * (1.0 - AX2 / (8.0 * NU + 16.0) * (1.0 - AX2 / (12.0 * NU + 36.0))));
       }else if((pow(L, 2.0) / AX) < 5.0 - 1){
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


void flag_spherbessel_basis(double *jell, const int ell, const double *nodes, int Nnodes)
{
	assert(Nnodes > 0);
	int i;
	for (i = 0; i < Nnodes; i++)
		jell[i] = j_ell(nodes[i], ell);
}




void flag_core_fourierbessel_analysis(complex double *flmn,
		const complex double *f,
		int L, double tau, int N)
{
	assert(L > 0);
	assert(N > 1);
	//const int alpha = ALPHA;
	int spin = 0;
	int verbosity = 0;
	int n;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
	int offset_lm, offset_r;

	complex double *flmr;
	flag_core_allocate_flmn(&flmr, L, N);

	for (n = 0; n < N; n++){
		//printf("> Analysis: layer %i on %i\n", n+1,N+1);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_forward_sov_conv_sym(flmr + offset_lm, f + offset_r, L, spin, dl_method, verbosity);
	}

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	assert(nodes != NULL);
	assert(weights != NULL);
	flag_spherlaguerre_sampling(nodes, weights, tau, N);
	//printf("> Mapped spherical Laguerre transform...");
	fflush(NULL);
	flag_spherlaguerre_mapped_analysis(flmn, flmr, nodes, weights, tau, N, flmsize);
	//printf("done\n");
	free(nodes);
	free(weights);
    free(flmr);

}

void flag_core_fourierbessel_synthesis(complex double *f,
		const complex double *flmn,
		const double *nodes, int Nnodes,
		int L, double tau, int N)
{
	assert(L > 0);
	assert(N > 1);
	assert(nodes != NULL);
	//const int alpha = ALPHA;
	int spin = 0;
	int verbosity = 0;
	int n, offset_lm, offset_r;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;

	complex double *flmr;
	flag_core_allocate_flmn(&flmr, L, Nnodes);
	//printf("> Mapped spherical Laguerre transform...");fflush(NULL);
	flag_spherlaguerre_mapped_synthesis(flmr, flmn, nodes, Nnodes, tau, N, flmsize);
	//printf("done\n");

	for (n = 0; n < Nnodes; n++){
		//printf("> Synthesis: layer %i on %i\n",n+1,Nnodes);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_inverse_sov_sym(f + offset_r, flmr + offset_lm, L, spin, dl_method, verbosity);
	}

    free(flmr);
}