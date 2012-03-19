// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * Compute forward FLAG transform (synthesis).
 *
 * Usage: 
 *   f = ...
 *     flag_synthesis_mex(flmn, L, N, reality);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int n, i, L, N, reality, flmn_m, flmn_n;
  double *flmn_real, *flmn_imag, *f_real, *f_imag;
  double *fr = NULL;
  complex double *flmn = NULL, *f = NULL;
  int ntheta, nphi, t, p;
  int iin = 0, iout = 0;

  /* Check number of arguments. */
  if(nrhs!=4) {
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:nrhs",
		      "Require four inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidOutput:nlhs",
		      "Require one output.");
  }

  /* Parse reality. */
  iin = 3;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:reality",
		      "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);

  /* Parse harmonic coefficients flm. */
  iin = 0;
  flmn_m = mxGetM(prhs[iin]);
  flmn_n = mxGetN(prhs[iin]);
  flmn_real = mxGetPr(prhs[iin]);
  flmn_imag = mxGetPi(prhs[iin]);
  flmn = (complex double*)malloc(flmn_m * flmn_n * sizeof(complex double));
  for (i=0; i<flmn_m*flmn_n; i++)
    flmn[i] = flmn_real[i];
  if( mxIsComplex(prhs[iin]) ) {
    flmn_imag = mxGetPi(prhs[iin]);
    for (i=0; i<flmn_m*flmn_n; i++)
      flmn[i] += I * flmn_imag[i];
  }

  /* Parse harmonic band-limit L. */
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit L must be positive integer.");
 /* Parse harmonic band-limit N. */
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit N must be positive integer.");

  if (flmn_m * flmn_n != L * L * N)
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:flmnSize",
		      "Invalid number of harmonic coefficients.");

  //printf("L = %i\n",L);
  //printf("N = %i\n",N);

  ntheta = L;
  nphi = 2 * L - 1;

  if (reality) {
	  fr = (double*)calloc(ntheta*nphi*N, sizeof(double));
	  flag_synthesis_real(fr, flmn, L, N);
  } else {
	  f = (complex double*)calloc(ntheta*nphi*N, sizeof(complex double));
	  flag_synthesis(f, flmn, L, N); 
  }

  //mexPrintf("flmn_m = %d; flmn_n = %d\n", flmn_m, flmn_n); 
  //mexPrintf("L = %d; N = %d; reality = %d\n", L, N, reality); 


  if (reality) {

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(N, ntheta*nphi, mxREAL);
    f_real = mxGetPr(plhs[iout]);
    for(n=0; n<N; n++)
      for(i=0; i<ntheta*nphi; i++) 
	        f_real[n*ntheta*nphi + i] = fr[n*ntheta*nphi + i];

  } else {

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(N, ntheta*nphi, mxCOMPLEX);
    f_real = mxGetPr(plhs[iout]);
    f_imag = mxGetPi(plhs[iout]);
    for(n=0; n<N; n++) {    
      for(i=0; i<ntheta*nphi; i++) {
      	  f_real[i*N + n] = creal(f[n*ntheta*nphi + i]);
	        f_imag[i*N + n] = cimag(f[n*ntheta*nphi + i]);
      }
    }

  }

  free(flmn);
  if (reality)
    free(fr);
  else
    free(f);


}