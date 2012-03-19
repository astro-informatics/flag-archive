// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#include <flag.h>
#include "mex.h"

/**
 * Compute forward FLAG transform (synthesis).
 *
 * Usage: 
 *   flmn = ...
 *     flag_analysis_mex(f, L, N, reality);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int i, L, N, n, f_m, f_n, reality;
  int f_is_complex;
  double *flmn_real = NULL, *flmn_imag = NULL;
  double *f_real = NULL, *f_imag = NULL;
  double *fr = NULL;
  complex double *flmn = NULL, *f = NULL;
  int ntheta, nphi, t, p;
  int iin = 0, iout = 0;

  /* Check number of arguments. */
  if(nrhs!=4) {
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:nrhs",
          "Require four inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidOutput:nlhs",
          "Require one output.");
  }

  /* Parse reality. */
  iin = 3;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:reality",
          "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);

  /* Parse function samples f. */
  iin = 0;
  if( !mxIsDouble(prhs[iin]) ) {
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:f",
          "Function values must be doubles.");
  }
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);  
  f_real = mxGetPr(prhs[iin]);  
  f_is_complex = mxIsComplex(prhs[iin]);
  f_imag = f_is_complex ? mxGetPi(prhs[iin]) : NULL;
  if (reality) {
    fr = (double*)malloc(f_m * f_n * sizeof(double));
    for(t=0; t<f_m; t++)
      for(p=0; p<f_n; p++)
        fr[t*f_n + p] = f_real[p*f_m + t];
  }
  else {
    f = (complex double*)malloc(f_m * f_n * sizeof(complex double));
    for(t=0; t<f_m; t++)
      for(p=0; p<f_n; p++)
        f[t*f_n + p] = f_real[p*f_m + t] 
          + I * (f_is_complex ? f_imag[p*f_m + t] : 0.0);
  }
  if (f_is_complex && reality)
    mexWarnMsgTxt("Running real transform but input appears to be complex (ignoring imaginary component).");
  if (!f_is_complex && !reality)
    mexWarnMsgTxt("Running complex transform on real signal (set reality flag to improve performance).");

  /* Parse harmonic band-limit L. */
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:bandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:bandLimit",
          "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit N must be positive integer.");


  //printf("L = %i\n",L);
  //printf("N = %i\n",N);

  ntheta = L;
  nphi = 2 * L - 1;
  flmn = (complex double*)calloc(L*L*N, sizeof(complex double));
  if (reality) {
    flag_analysis_real(flmn, fr, L, N);
  } else {
    flag_analysis(flmn, f, L, N); 
  }

  //mexPrintf("f_m = %d; f_n = %d\n", f_m, f_n); 
  //mexPrintf("L = %d; N = %d; reality = %d\n", L, N, reality); 


  /* Copy flm to output. */
  plhs[iout] = mxCreateDoubleMatrix(N, L*L, mxCOMPLEX);
  flmn_real = mxGetPr(plhs[iout]);
  flmn_imag = mxGetPi(plhs[iout]);
  for (i=0; i<L*L*N; i++) {
    flmn_real[i] = creal(flmn[i]);
    flmn_imag[i] = cimag(flmn[i]);;
  }

  /* Free memory. */
  free(flmn);
  if (reality)
    free(fr);
  else
    free(f);
}