// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * MATLAB interface: flag_get_tau_mex.
 * This function for internal use only.
 * Compute scaling factor for SLAG transform.
 *
 * Usage: 
 *   tau = flag_get_tau_mex(N, R);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int N;
  double R;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=2) {
    mexErrMsgIdAndTxt("flag_get_tau_mex:InvalidInput:nrhs",
          "Require two inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("flag_get_tau_mex:InvalidOutput:nlhs",
          "Require one output.");
  }

  // Parse harmonic band-limit N
  iin = 0;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_get_tau_mex:InvalidInput:bandLimit",
          "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("flag_get_tau_mex:InvalidInput:bandLimitNonInt",
          "Harmonic band-limit N must be positive integer.");
  // Parse radial limit R
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_get_tau_mex:InvalidInput:radialLimit",
          "Radial limit R must be positive double.");
  }
  R = mxGetScalar(prhs[iin]);
  if (R <= 0)
    mexErrMsgIdAndTxt("flag_get_tau_mex:InvalidInput:radialLimit",
          "Radial limit R must be positive double.");

  double tau = flag_spherlaguerre_tau(R, N);

  // Copy result to matlab output
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, 1,mxREAL);
  double *tau_out;
  tau_out = mxGetPr(plhs[iout]);
  tau_out[0] = tau;

}
