// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * MATLAB interface: slag_sampling.
 * This function for internal use only.
 * Compute Spherical Laguerre sampling
 *
 * Usage: 
 *   nodes = ...
 *     slag_sampling_mex(N, R);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int n, N;
  double R;
  double *nodes_out, *weights;
  double *nodes = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=2) {
    mexErrMsgIdAndTxt("slag_sampling_mex:InvalidInput:nrhs",
          "Require two inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("slag_sampling_mex:InvalidOutput:nlhs",
          "Require one output.");
  }

 // Parse harmonic band-limit N
  iin = 0;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_sampling_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("slag_sampling_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit N must be positive integer.");

  // Parse harmonic band-limit R
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_sampling_mex:InvalidInput:Rlimit",
          "Radial limit R must be positive real.");
  }
  R = mxGetScalar(prhs[iin]);
  if ( R <= 0.0 )
    mexErrMsgIdAndTxt("slag_sampling_mex:InvalidInput:RLimitNonInt",
          "Radial limit R must be positive real.");

  nodes = (double*)calloc(N, sizeof(double));

  // Compute nodes
  weights = (double*)calloc(N, sizeof(double));
  flag_spherlaguerre_sampling(nodes, weights, R, N);
  free(weights);

  //printf("N = %i\n",N);
  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, N, mxREAL);
  nodes_out = mxGetPr(plhs[iout]);
  for(n=0; n<N; n++)
   nodes_out[n] = nodes[n];

  free(nodes);

}
