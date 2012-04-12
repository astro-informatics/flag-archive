// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAG_SPHERBESSEL
#define FLAG_SPHERBESSEL

void flag_spherlaguerre2spherbessel(double *flk, const double *fn, double *kvalues, int Nk, int N, int ell, double tau);

void flag_spherbessel_approx(double *flk, const double *f, double *kvalues, int Nk, double *nodes, int Nnodes, int ell);

void flag_fourierlaguerre2fourierbessel(complex double *flmk, complex double *flmn, double *kvalues, int Nk, int N, int L, double tau);

void flag_sbesselslag(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau);

double flag_mujlk(int j, int ell, double k, double tau);

double j_ell(double X, int l);

#endif