// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * MATLAB interface: spherbessel_basis.
 * This function for internal use only.
 *
 * Usage: 
 *   j_ell(radii) = spherbessel_basis_mex(ell, nodes);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int n, i, ell, nodes_m, nodes_n;
  double *f_real, *nodes_in;
  double *jell = NULL, *nodes = NULL;
  int iin = 0, iout = 0;

  // Check number of arguments
  if(nrhs!=2) {
    mexErrMsgIdAndTxt("spherbessel_basis_mex:InvalidInput:nrhs",
          "Require three inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("spherbessel_basis_mex:InvalidOutput:nlhs",
          "Require one output.");
  }

 // Parse spherical Bessel mode ell
  iin = 0;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("spherbessel_basis_mex:InvalidInput:Nmode",
		      "spherical Laguerre N must be integer.");
  }
  ell = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)ell || ell < 0)
    mexErrMsgIdAndTxt("spherbessel_basis_mex:InvalidInput:NmodeNonInt",
		      "spherical Bessel ell must be positive integer.");

  // Parse nodes
  iin = 1;
  nodes_m = mxGetM(prhs[iin]);
  nodes_n = mxGetN(prhs[iin]);
  int Nnodes = nodes_m*nodes_n;
  nodes_in = mxGetPr(prhs[iin]);
  nodes = (double*)calloc(Nnodes, sizeof(double));
  for (i=0; i<Nnodes; i++)
    nodes[i] = nodes_in[i];

  jell = (double*)calloc(Nnodes, sizeof(double));
  flag_spherbessel_basis(jell, ell, nodes, Nnodes);

  iout = 0;
  plhs[iout] = mxCreateDoubleMatrix(nodes_m, nodes_n, mxREAL);
  f_real = mxGetPr(plhs[iout]);
  for(n=0; n<nodes_m*nodes_n; n++)
    f_real[n] = jell[n];
    
  free(jell);
  free(nodes);

}
