// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * MATLAB interface: slag_basis.
 * This function for internal use only.
 *
 * Usage: 
 *   K_N(radii) = slag_basis_mex(N, nodes, tau);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int n, i, N, Nnodes, nodes_m, nodes_n;
  double tau;
  double *f_real, *nodes_in;
  double *KN = NULL, *nodes = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=3) {
    mexErrMsgIdAndTxt("slag_basis_mex:InvalidInput:nrhs",
          "Require three inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("slag_basis_mex:InvalidOutput:nlhs",
          "Require one output.");
  }

 // Parse spherical Laguerre mode N
  iin = 0;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_basis_mex:InvalidInput:Nmode",
		      "spherical Laguerre N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N < 0)
    mexErrMsgIdAndTxt("slag_basis_mex:InvalidInput:NmodeNonInt",
		      "spherical Laguerre N must be positive integer.");

  // Parse nodes
  iin = 1;
  nodes_m = mxGetM(prhs[iin]);
  nodes_n = mxGetN(prhs[iin]);
  Nnodes = nodes_m*nodes_n;
  nodes_in = mxGetPr(prhs[iin]);
  nodes = (double*)calloc(Nnodes, sizeof(double));
  for (i=0; i<Nnodes; i++)
    nodes[i] = nodes_in[i];

  // Parse harmonic band-limit tau
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("slag_basis_mex:InvalidInput:Rlimit",
          "Scale factor tau must be positive real.");
  }
  tau = mxGetScalar(prhs[iin]);
  if ( tau <= 0 )
    mexErrMsgIdAndTxt("slag_basis_mex:InvalidInput:RLimitNonInt",
          "Scale factor tau must be positive real.");
  
  KN = (double*)calloc(Nnodes*N, sizeof(double));
  flag_spherlaguerre_basis(KN, N, nodes, Nnodes, tau);

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(Nnodes, N, mxREAL);
  f_real = mxGetPr(plhs[iout]);
  for(n=0; n<N; n++)
    for(i=0; i<Nnodes; i++)
      f_real[n*Nnodes+i] = KN[i*N+n];
    
  free(KN);
  free(nodes);

}
