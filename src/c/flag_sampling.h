
#ifndef FLAG_SAMPLING
#define FLAG_SAMPLING

enum ssht_methods { 
	MW, 
	MWSS, 
	GL, 
	DH 
};

void flag_allocate_sampling(double **rs, double **thetas, double **phis, double **laguweights, double R, int L, int N, enum ssht_methods method);

void flag_deallocate_sampling(double *rs, double *thetas, double *phis, double *laguweights);

void ssht_allocate_sampling(double **thetas, double **phis, int L, enum ssht_methods method);

void ssht_sampling(double *thetas, double *phis, int L, enum ssht_methods method);

void flag_sampling(double *rs, double *thetas, double *phis, double *laguweights, double R, int L, int N, enum ssht_methods method);

#endif