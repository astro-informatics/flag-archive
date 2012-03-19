// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAG_SPHERBESSEL
#define FLAG_SPHERBESSEL

void flag_spherlaguerre2spherbessel(complex double *flk, complex double *fn, double *kvalues, int Nk, int N, int ell, double tau);

double flag_mujlk(int j, int l, double k, double tau);

#endif