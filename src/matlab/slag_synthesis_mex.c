// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * MATLAB interface: slag_synthesis.
 * This function for internal use only.
 * Compute forward 1D SLAG transform (synthesis).
 *
 * Usage: 
 *   [f, nodes] = ...
 *     slag_synthesis_mex(fn, N, R, nodes);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int n, i, N, fn_m, fn_n, nodes_m, nodes_n;
  double R;
  double *fn_real, *f_real, *nodes_in, *nodes_out, *weights;
  double *fn = NULL, *f = NULL, *nodes = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=4) {
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:nrhs",
          "Require four inputs.");
  }
  if(nlhs!=2) {
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidOutput:nlhs",
          "Require two outputs.");
  }

  // Parse harmonic coefficients fn
  iin = 0;
  fn_m = mxGetM(prhs[iin]);
  fn_n = mxGetN(prhs[iin]);
  fn_real = mxGetPr(prhs[iin]);
  fn = (double*)malloc(fn_m * fn_n * sizeof(double));
  for (i=0; i<fn_m*fn_n; i++)
    fn[i] = fn_real[i];
  

 // Parse harmonic band-limit N
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit N must be positive integer.");

  if (fn_m * fn_n != N)
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:fnSize",
		      "Invalid number of harmonic coefficients.");

  // Parse harmonic band-limit R
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:Rlimit",
          "Radial limit R must be positive real.");
  }
  R = mxGetScalar(prhs[iin]);
  if ( R <= 0 )
    mexErrMsgIdAndTxt("slag_synthesis_mex:InvalidInput:RLimitNonInt",
          "Radial limit R must be positive real.");
  

  // Parse nodes
  iin = 3;
  nodes_m = mxGetM(prhs[iin]);
  nodes_n = mxGetN(prhs[iin]);
  if( nodes_m*nodes_n > 1 ){

    f = (double*)calloc(nodes_m*nodes_n, sizeof(double));

    nodes_in = mxGetPr(prhs[iin]);
    nodes = (double*)calloc(nodes_m*nodes_n, sizeof(double));
    for (i=0; i<nodes_m*nodes_n; i++)
      nodes[i] = nodes_in[i];

    // Run spherical Laguerre synthesis
    flag_spherlaguerre_synthesis(f, fn, nodes, nodes_m*nodes_n, N);

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(nodes_m, nodes_n, mxREAL);
    f_real = mxGetPr(plhs[iout]);
    for(n=0; n<nodes_m*nodes_n; n++)
      f_real[n] = f[n];
    iout = 1;
    plhs[iout] = mxCreateDoubleMatrix(nodes_m, nodes_n, mxREAL);
    nodes_out = mxGetPr(plhs[iout]);
    for(n=0; n<nodes_n*nodes_m; n++)
      nodes_out[n] = nodes[n];

  }else{

    f = (double*)calloc(N+1, sizeof(double));

    weights = (double*)calloc(N+1, sizeof(double));
    nodes = (double*)calloc(N+1, sizeof(double));
    flag_spherlaguerre_sampling(nodes, weights, R, N);
    free(weights);

    // Run spherical Laguerre synthesis
    flag_spherlaguerre_synthesis(f, fn, nodes, N+1, N);

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(1, N+1, mxREAL);
    f_real = mxGetPr(plhs[iout]);
    for(n=0; n<N+1; n++)
      f_real[n] = f[n];
    iout = 1;
    plhs[iout] = mxCreateDoubleMatrix(1, N+1, mxREAL);
    nodes_out = mxGetPr(plhs[iout]);
    for(n=0; n<N+1; n++)
      nodes_out[n] = nodes[n];

  }

  free(fn);
  free(nodes);
  free(f);

}
