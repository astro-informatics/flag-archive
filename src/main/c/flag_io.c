// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flag.h"
#include <math.h>
#include <stdlib.h>
#include <complex.h> 
#include <fftw3.h> 
#include <ssht.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

void flag_write_f_real(const double *f_real, int L, int N, char *outfile)
{
  assert(L > 0);
  assert(N > 1);
  int l, n, frsize = ssht_fr_size(L);

  FILE *file; 
  file = fopen(outfile,"a+");

  //fprintf(file, "%i", L);
  //fprintf(file, "%i", N);

  for (n = 0; n < N; n++){
    for(l=0; l<frsize; l++){
      fprintf(file, "%f ", f_real[l+n*frsize]);
    }
    fprintf(file, "\n");
  }

  fclose(file);
}

void flag_write_f(const complex double *f, int L, int N, char *outfile)
{
  assert(L > 0);
  assert(N > 1);
  int l, n, frsize = ssht_fr_size(L);

  FILE *file; 
  file = fopen(outfile,"a+");

  //fprintf(file, "%i", L);
  //fprintf(file, "%i", N);

  for (n = 0; n < N; n++){
    for(l=0; l<frsize; l++){
      fprintf(file, "%f %f ", creal(f[l+n*frsize]), cimag(f[l+n*frsize]));
    }
    fprintf(file, "\n");
  }

  fclose(file);
}
