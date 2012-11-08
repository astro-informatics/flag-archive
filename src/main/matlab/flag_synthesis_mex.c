// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <flag.h>
#include "mex.h"

/**
 * MATLAB interface: flag_synthesis.
 * This function for internal use only.
 * Compute forward FLAG transform (synthesis).
 *
 * Usage: 
 *   f = ...
 *     flag_synthesis_mex(flmn, L, N, nodes, R, reality);
 *
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int f_is_complex, n, i, L, N, reality, flmn_m, flmn_n, Nnodes;
  double *flmn_real, *flmn_imag, *f_real, *f_imag;
  double *fr = NULL;
  complex double *flmn = NULL, *f = NULL;
  int ntheta, nphi;
  int iin = 0, iout = 0;


  // Check number of arguments
  if(nrhs!=6) {
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:nrhs",
		      "Require six inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidOutput:nlhs",
		      "Require one output.");
  }

  // Parse reality flag
  iin = 5;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:reality",
		      "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);


  // Parse harmonic coefficients flm
  iin = 0;
  flmn_m = mxGetM(prhs[iin]);
  flmn_n = mxGetN(prhs[iin]);
  flmn_real = mxGetPr(prhs[iin]);
  flmn_imag = mxGetPi(prhs[iin]);
  f_is_complex = mxIsComplex(prhs[iin]);
  flmn = (complex double*)malloc(flmn_m * flmn_n * sizeof(complex double));
  for(n=0; n<flmn_m; n++)  
    for(i=0; i<flmn_n; i++) 
      flmn[n*flmn_n+i] = flmn_real[i*flmn_m+n];
  for(n=0; n<flmn_m; n++)  
    for(i=0; i<flmn_n; i++) 
      flmn[n*flmn_n+i] += I * (f_is_complex ? flmn_imag[i*flmn_m+n] : 0.0);


  // Parse harmonic band-limit L
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit L must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit L must be positive integer.");
 /* Parse harmonic band-limit N. */
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:LbandLimit",
		      "Harmonic band-limit N must be integer.");
  }
  N = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)N || N <= 0)
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit N must be positive integer.");

  if (flmn_m * flmn_n != L * L * N)
    mexErrMsgIdAndTxt("flag_synthesis_mex:InvalidInput:flmnSize",
		      "Invalid number of harmonic coefficients.");

  ntheta = L;
  nphi = 2 * L - 1;

  // Parse harmonic band-limit R
  iin = 4;
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

  // Parse nodes
  int nodes_m, nodes_n;
  iin = 3;
  nodes_m = mxGetM(prhs[iin]);
  nodes_n = mxGetN(prhs[iin]);
  double *nodes_real, *nodes, *weights;


  if(nodes_m*nodes_n > 1){
    nodes_real = mxGetPr(prhs[iin]);
    nodes = (double*)malloc(nodes_m*nodes_n * sizeof(double));  
    for(n=0; n<nodes_m*nodes_n; n++)  
      nodes[n] = nodes_real[n];
    Nnodes = nodes_m*nodes_n;
  }else{
    Nnodes = N ;
    nodes = (double*)calloc(Nnodes, sizeof(double));
    weights  = (double*)calloc(Nnodes, sizeof(double));
    flag_spherlaguerre_sampling(nodes, weights, R, N);
    free(weights);
  }


  if (reality) {
  	fr = (double*)calloc(ntheta*nphi*Nnodes, sizeof(double));
  	flag_core_synthesis_real(fr, flmn, nodes, Nnodes, L, N);
  } else {
  	 f = (complex double*)calloc(ntheta*nphi*Nnodes, sizeof(complex double));
  	 flag_core_synthesis(f, flmn, nodes, Nnodes, L, N); 
  }


  if (reality) {

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(Nnodes, ntheta*nphi, mxREAL);
    f_real = mxGetPr(plhs[iout]);
    for(n=0; n<Nnodes; n++)
      for(i=0; i<ntheta*nphi; i++) 
	        f_real[i*Nnodes + n] = fr[n*ntheta*nphi + i];

  } else {

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(Nnodes, ntheta*nphi, mxCOMPLEX);
    f_real = mxGetPr(plhs[iout]);
    f_imag = mxGetPi(plhs[iout]);
    for(n=0; n<Nnodes; n++) {    
      for(i=0; i<ntheta*nphi; i++) {
      	  f_real[i*Nnodes + n] = creal(f[n*ntheta*nphi + i]);
	        f_imag[i*Nnodes + n] = cimag(f[n*ntheta*nphi + i]);
      }
    }

  }

  free(flmn);
  if (reality)
    free(fr);
  else
    free(f);

  free(nodes);

}
