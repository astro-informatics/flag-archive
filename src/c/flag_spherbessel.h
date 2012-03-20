// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAG_SPHERBESSEL
#define FLAG_SPHERBESSEL

void flag_spherlaguerre2spherbessel(double *flk, double *fn, double *kvalues, int Nk, int N, int ell, double tau);

void flag_sbesselslag(double *sbesselslag, int ell, double *kvalues, int Nk, int N, double tau);

double flag_mujlk(int j, int ell, double k, double tau);

#endif