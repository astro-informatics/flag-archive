// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * MATLAB interface: flag_analysis.
 * This function for internal use only.
 * Compute forward FLAG transform (synthesis).
 *
 * Usage: 
 *   flmn = ...
 *     flag_analysis_mex(f, L, N, R, reality);
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
  if(nrhs!=5) {
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:nrhs",
          "Require five inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidOutput:nlhs",
          "Require one output.");
  }

  // Parse reality flag
  iin = 4;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("flag_analysis_mex:InvalidInput:reality",
          "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);

  // Parse function samples f
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


  // Parse harmonic band-limit L
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

  // Parse harmonic band-limit R
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:Rlimit",
          "Radial limit R must be positive real.");
  }
  double R = mxGetScalar(prhs[iin]);
  if ( R <= 0 )
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:RLimitNonInt",
          "Radial limit R must be positive real.");


  ntheta = L;
  nphi = 2 * L - 1;
  flmn = (complex double*)calloc(L*L*N, sizeof(complex double));
  if (reality) {
    flag_core_analysis_real(flmn, fr, R, L, N);
  } else {
    flag_core_analysis(flmn, f, R, L, N); 
  }

  // Copy flm to output
  plhs[iout] = mxCreateDoubleMatrix(N, L*L, mxCOMPLEX);
  flmn_real = mxGetPr(plhs[iout]);
  flmn_imag = mxGetPi(plhs[iout]);
  for(n=0; n<N; n++) {    
    for(i=0; i<L*L; i++) {
      flmn_real[i*N+n] = creal(flmn[n*L*L+i]);
      flmn_imag[i*N+n] = cimag(flmn[n*L*L+i]);;
    }
  }

  free(flmn);
  if (reality)
    free(fr);
  else
    free(f);
}
