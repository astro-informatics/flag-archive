// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAG_SAMPLING
#define FLAG_SAMPLING

/* FOR FUTURE IMPROVEMENTS // multi-scheme support
enum ssht_methods { 
	MW, 
	MWSS, 
	GL, 
	DH 
};
*/

void flag_allocate_sampling(double **rs, double **thetas, double **phis, double **laguweights, double R, int L, int N);

void ssht_allocate_sampling(double **thetas, double **phis, int L);

void ssht_sampling(double *thetas, double *phis, int L);

void flag_sampling(double *rs, double *thetas, double *phis, double *laguweights, double R, int L, int N);

#endif