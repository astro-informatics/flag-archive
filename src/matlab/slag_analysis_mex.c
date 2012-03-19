// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * Compute backward 1D SLAG transform (analysis).
 *
 * Usage: 
 *   fn = ...
 *     slag_analysis_mex(f, N, R);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int n, i, N, f_m, f_n;
  double R;
  double *fn_real, *f_real;
  double *fn = NULL, *f = NULL, *nodes = NULL, *weights = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=3) {
    mexErrMsgIdAndTxt("slag_analysis_mex:InvalidInput:nrhs",
          "Require three inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("slag_analysis_mex:InvalidOutput:nlhs",
          "Require one output.");
  }

  // Parse harmonic coefficients fn
  iin = 0;
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);
  f_real = mxGetPr(prhs[iin]);
  f = (double*)malloc(f_m * f_n * sizeof(double));
  for (i=0; i<f_m*f_n; i++)
    f[i] = f_real[i];
  

 // Parse harmonic band-limit N
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_analysis_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("slag_analysis_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit N must be positive integer.");

  if (f_m * f_n != N)
    mexErrMsgIdAndTxt("slag_analysis_mex:InvalidInput:fnSize",
		      "Invalid number of harmonic coefficients.");

  // Parse harmonic band-limit R
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_analysis_mex:InvalidInput:Rlimit",
          "Radial limit R must be positive real.");
  }
  R = mxGetScalar(prhs[iin]);
  if ( R <= 0)
    mexErrMsgIdAndTxt("slag_analysis_mex:InvalidInput:RLimitNonInt",
          "Radial limit R must be positive real.");

  nodes = (double*)calloc(N, sizeof(double));

  // Parse nodes
  weights = (double*)calloc(N, sizeof(double));
  flag_spherlaguerre_sampling(nodes, weights, R, N);

  // Run spherical Laguerre analysis
  fn = (double*)calloc(N, sizeof(double));
	flag_spherlaguerre_analysis(fn, f, nodes, weights, N);

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, N, mxREAL);
  fn_real = mxGetPr(plhs[iout]);
  for(n=0; n<N; n++)
	  fn_real[n] = fn[n];

  free(fn);
  free(nodes);
  free(f);
  free(weights);

}
