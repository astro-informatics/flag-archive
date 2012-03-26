
#ifndef FLAG_SPHERLAGUERRE
#define FLAG_SPHERLAGUERRE

void flag_spherlaguerre_quadrature(double *roots, double *weights, int N);

double flag_spherlaguerre_tau(double R, int N);

void flag_spherlaguerre_sampling(double *nodes, double *weights, double R, int N);

void flag_spherlaguerre_analysis(double *fn, const double *f, const double *nodes, const double *weights, int N);

void flag_spherlaguerre_synthesis(double *f, const double *fn, const double *nodes, int Nnodes, int N);

void flag_allocate_spherlaguerre_sampling(double **nodes, double **weights, int N);

void flag_deallocate_spherlaguerre_sampling(double *nodes, double *weights);

void flag_mapped_spherlaguerre_analysis(complex double *fn, const complex double *f, const double *weights, const double *nodes, int N, int mapsize);

void flag_mapped_spherlaguerre_synthesis(complex double *f, const complex double *fn, const double *nodes, int Nnodes, int N, int mapsize);

#endif
