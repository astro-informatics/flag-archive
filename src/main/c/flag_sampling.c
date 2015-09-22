// FLAG package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "flag.h"
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <ssht.h>
#include <assert.h>

void allocate_ssht_sampling(double **thetas, double **phis, int L)
{
  assert(L > 0);
	int nphi, ntheta;

  ntheta = ssht_sampling_mw_ntheta(L);
  nphi = ssht_sampling_mw_nphi(L);

  	/* FOR FUTURE IMPROVEMENTS // multi-scheme support
  	switch (method)
  	{
  		case MW:
  			ntheta = ssht_sampling_mw_ntheta(L);
    		nphi = ssht_sampling_mw_nphi(L);
    		break;
	  	case MWSS:
	  		ntheta = ssht_sampling_mw_ss_ntheta(L);
	    	nphi = ssht_sampling_mw_ss_nphi(L);
	    	break;
	    case GL:
	    	ntheta = ssht_sampling_gl_ntheta(L);
    		nphi = ssht_sampling_gl_nphi(L);
    		break;
      	case DH:
      		ntheta = ssht_sampling_dh_ntheta(L);
    		nphi = ssht_sampling_dh_nphi(L);
    		break;
	  }
	  */

	  *thetas = (double*)calloc(ntheta, sizeof(double));
	  *phis = (double*)calloc(nphi, sizeof(double));

}

void flag_sampling_allocate(double **rs, double **thetas, double **phis, double **laguweights, int L, int N)
{
  assert(L > 0);
  assert(N > 1);
  allocate_ssht_sampling(thetas, phis, L);
  flag_spherlaguerre_allocate_sampling(rs, laguweights, N);
}

void ssht_sampling_mw(double *thetas, double *phis, int L)
{
  assert(L > 0);
	int t, p, nphi, ntheta;

	ntheta = ssht_sampling_mw_ntheta(L);
  nphi = ssht_sampling_mw_nphi(L);
  for (t=0; t<ntheta; t++)
	 	thetas[t] = ssht_sampling_mw_t2theta(t, L);
	for (p=0; p<nphi; p++)
	  phis[p] = ssht_sampling_mw_p2phi(p, L);
}

void flag_sampling(double *rs, double *thetas, double *phis, double *laguweights, double tau, int L, int N)
{
  assert(L > 0);
  assert(N > 1);
  assert(tau > 0.0);
	flag_spherlaguerre_sampling(rs, laguweights, tau, N);
	ssht_sampling_mw(thetas, phis, L);
}
