// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAG_CORE
#define FLAG_CORE

void flag_random_f(complex double *f, int L, int N, int seed);

void flag_random_flmn(complex double *flmn, int L, int N, int seed);

void flag_allocate_flmn(complex double **flmn, int L, int N);

int ssht_fr_size(int L);

int ssht_flm_size(int L);

void flag_analysis(complex double *flmn, const complex double *f, int L, int N);

void flag_synthesis(complex double *f, const complex double *flmn, int L, int N);

void flag_analysis_real(complex double *flmn, const double *f, int L, int N);

void flag_synthesis_real(double *f, const complex double *flmn, int L, int N);

#endif