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

void flag_allocate_flmn(complex double **flmn, int L, int N)
{
	assert(L > 0);
	assert(N > 1);
	int flmsize = ssht_flm_size(L);
	long totalsize = N*flmsize;
	*flmn = (complex double*)calloc(totalsize, sizeof(complex double));
	assert(flmn != NULL);
}

void flag_allocate_f(complex double **f, int L, int N)
{
	assert(L > 0);
	assert(N > 1);
	int frsize = ssht_fr_size_mw(L);
	long totalsize = (N) * frsize;
	*f = (complex double*)calloc(totalsize, sizeof(complex double));
	assert(f != NULL);
}

void flag_allocate_f_real(double **f, int L, int N)
{
	assert(L > 0);
	assert(N > 1);
	int frsize = ssht_fr_size_mw(L);
	long totalsize = (N) * frsize;
	*f = (double*)calloc(totalsize, sizeof(double));
	assert(f != NULL);
}

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

int flag_flmn_size(int L, int N)
{// In case we want to extend to various sampling schemes
	assert(L > 0);
	assert(N > 1);
	return ssht_flm_size(L)*N;
}

int flag_f_size_mw(int L, int N)
{// In case we want to extend to various sampling schemes
	assert(L > 0);
	assert(N > 1);
	return ssht_fr_size_mw(L) * (N);
}

void flag_analysis(complex double *flmn, 
		const complex double *f, 
		double R, int L, int N)
{
	assert(L > 0);
	assert(N > 1);
	//const int alpha = 2;
	int spin = 0;
	int verbosity = 0;
	int n;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
	int offset_lm, offset_r;
	
	complex double *flmr;
	flag_allocate_flmn(&flmr, L, N);
	
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
	flag_spherlaguerre_sampling(nodes, weights, R, N);
	//printf("> Mapped spherical Laguerre transform...");
	fflush(NULL);
	flag_mapped_spherlaguerre_analysis(flmn, flmr, nodes, weights, N, flmsize);
	//printf("done\n");
	free(nodes);
	free(weights);
    free(flmr);

}

void flag_synthesis(complex double *f, 
		const complex double *flmn, 
		const double *nodes, int Nnodes,
		int L, int N)
{
	assert(L > 0);
	assert(N > 1);
	assert(nodes != NULL);
	//const int alpha = 2;
	int spin = 0;
	int verbosity = 0;
	int n, offset_lm, offset_r;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;

	complex double *flmr;
	flag_allocate_flmn(&flmr, L, Nnodes);
	//printf("> Mapped spherical Laguerre transform...");fflush(NULL);
	flag_mapped_spherlaguerre_synthesis(flmr, flmn, nodes, Nnodes, N, flmsize);
	//printf("done\n");
	
	for (n = 0; n < Nnodes; n++){
		//printf("> Synthesis: layer %i on %i\n",n+1,Nnodes);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_inverse_sov_sym(f + offset_r, flmr + offset_lm, L, spin, dl_method, verbosity);
	}

    free(flmr);
}

void flag_analysis_real(complex double *flmn, 
		const double *f, double R, int L, int N)
{
	assert(L > 0);
	assert(N > 1);
	int verbosity = 0;
	int n, offset_lm, offset_r;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;

	complex double *flmr;
	flag_allocate_flmn(&flmr, L, N);

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
	flag_spherlaguerre_sampling(nodes, weights, R, N);
	//printf("> Mapped spherical Laguerre transform...");
	fflush(NULL);
	flag_mapped_spherlaguerre_analysis(flmn, flmr, nodes, weights, N, flmsize);
	//printf("done\n");
	free(nodes);
	free(weights);
    free(flmr);
}

void flag_synthesis_real(double *f, 
		const complex double *flmn, 
		const double *nodes, int Nnodes, 
		int L, int N)
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
	flag_allocate_flmn(&flmr, L, Nnodes);
	flag_mapped_spherlaguerre_synthesis(flmr, flmn, nodes, Nnodes, N, flmsize);
	//printf("done\n");

	for (n = 0; n < Nnodes; n++){
		//printf("> Synthesis: layer %i on %i\n",n+1,N);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_inverse_sov_sym_real(f + offset_r, flmr + offset_lm, L, dl_method, verbosity);
	}

    free(flmr);
}