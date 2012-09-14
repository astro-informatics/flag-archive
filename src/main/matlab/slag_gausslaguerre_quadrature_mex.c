// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * MATLAB interface: slag_gausslaguerre_quadrature.
 * This function for internal use only.
 *
 * Usage: 
 *   [nodes, weights] = ...
 *     slag_gausslaguerre_quadrature_mex(N, alpha);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int n, N, alpha, nodes_m, nodes_n;
  double *nodes_out, *weights_out;
  double *weights = NULL, *nodes = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=2) {
    mexErrMsgIdAndTxt("slag_gausslaguerre_quadrature_mex:InvalidInput:nrhs",
          "Require two inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("slag_gausslaguerre_quadrature_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

 // Parse band-limit N
  iin = 0;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_gausslaguerre_quadrature_mex:InvalidInput:",
		      "Band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("slag_gausslaguerre_quadrature_mex:InvalidInput:",
		      "Band-limit N must be positive integer.");

  // Parse order alpha
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_gausslaguerre_quadrature_mex:InvalidInput:",
          "Order alpha must be integer.");
  }
  alpha = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)alpha || alpha < 0)
    mexErrMsgIdAndTxt("slag_gausslaguerre_quadrature_mex:InvalidInput:",
          "Order alpha must be positive integer.");


  nodes = (double*)calloc(N, sizeof(double));
  weights = (double*)calloc(N, sizeof(double));

  // Compute quadrature
  flag_spherlaguerre_quadrature(nodes, weights, N, alpha);

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(1, N, mxREAL);
  nodes_out = mxGetPr(plhs[iout]);
  for(n=0; n<N; n++)
    nodes_out[n] = nodes[n];
  
  iout = 1;
  plhs[iout] = mxCreateDoubleMatrix(1, N, mxREAL);
  weights_out = mxGetPr(plhs[iout]);
  for(n=0; n<N; n++)
    weights_out[n] = weights[n];

  free(nodes);
  free(weights);

}
