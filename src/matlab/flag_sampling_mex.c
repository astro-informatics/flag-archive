// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * MATLAB interface: flag_sampling.
 * This function for internal use only.
 * Compute FLAG sampling scheme.
 *
 * Usage: 
 *   [rs, thetas, phis] = ...
 *     flag_sampling(L, N, R);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int L, N, n;
  double R;
  double *rs = NULL, *tempweights = NULL, *thetas = NULL, *phis = NULL;
  double *rs_out, *thetas_out, *phis_out;
  int nthetas, nphis, nrs;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=3) {
    mexErrMsgIdAndTxt("flag_sampling_mex:InvalidInput:nrhs",
          "Require three inputs.");
  }
  if(nlhs!=3) {
    mexErrMsgIdAndTxt("flag_sampling_mex:InvalidOutput:nlhs",
          "Require three outputs.");
  }

  // Parse harmonic band-limit L
  iin = 0;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_sampling_mex:InvalidInput:bandLimit",
          "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("flag_sampling_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit L must be positive integer.");

  // Parse harmonic band-limit N
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_sampling_mex:InvalidInput:bandLimit",
          "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("flag_sampling_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit N must be positive integer.");
  // Parse radial limit R
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_sampling_mex:InvalidInput:radialLimit",
          "Radial limit R must be positive double.");
  }
  R = mxGetScalar(prhs[iin]);
  if (R <= 0)
    mexErrMsgIdAndTxt("flag_sampling_mex:InvalidInput:radialLimit",
          "Radial limit R must be positive double.");

  // MW sampling
  nrs = N;
  nthetas = L;
  nphis = 2 * L - 1;

  // Allocate arrays
  rs = (double*)calloc(nrs, sizeof(double));
  tempweights = (double*)calloc(nrs, sizeof(double));
  thetas = (double*)calloc(nthetas, sizeof(double));
  phis = (double*)calloc(nphis, sizeof(double));


  // Execute flag_sampling function
  flag_sampling(rs, thetas, phis, tempweights, R, L, N);
  

  // Copy result to matlab output
  // Rs
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, nrs,mxREAL);
  rs_out = mxGetPr(plhs[iout]);
  for (n=0; n<nrs; n++)
    rs_out[n] = rs[n];
  // Thetas
  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(1, nthetas, mxREAL);
  thetas_out = mxGetPr(plhs[iout]);
  for (n=0; n<nthetas; n++)
    thetas_out[n] = thetas[n];
  // Phis
  iout = 2;
  plhs[iout] = mxCreateDoubleMatrix(1, nphis, mxREAL);
  phis_out = mxGetPr(plhs[iout]);
  for (n=0; n<nphis; n++)
    phis_out[n] = phis[n];

  // Free memory
  free(rs);
  free(tempweights);
  free(thetas);
  free(phis);

}
