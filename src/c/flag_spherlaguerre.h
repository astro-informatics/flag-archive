
#ifndef FLAG_SPHERLAGUERRE
#define FLAG_SPHERLAGUERRE

void flag_spherlaguerre_quadrature(double *roots, double *weights, int N);

double flag_spherlaguerre_tau(double R, int N);

void flag_spherlaguerre_sampling(double *nodes, double *weights, double R, int N);

void flag_spherlaguerre_analysis(double *fn, const double *f, const double *nodes, const double *weights, int N);

void flag_spherlaguerre_synthesis(double *f, const double *fn, const double *nodes, int N);

void flag_spherlaguerre_sampling_allocate(double **nodes, double **weights, int N);

void flag_spherlaguerre_sampling_deallocate(double *nodes, double *weights);

#endif
