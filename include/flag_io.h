// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAG_IO
#define FLAG_IO

#include <complex.h> 

void flag_write_f(const complex double *f, int L, int N, char *outfile);
void flag_write_f_real(const double *f_real, int L, int N, char *outfile);

void flag_read_f(const complex double *f, int L, int N, char *infile);
void flag_read_f_real(const double *f_real, int L, int N, char *infile);

#endif