// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ssht.h>

#include "flag_sampling.h"
#include "flag_spherlaguerre.h"

/*!  
 * Compute the exact Fourier-Laguerre sampling.
 *
 * \param[out] f Sampling on the 3-ball.
 * \param[in] L Harmonic band-limit.
 * \param[in] N Laguerre band-limit.
 * \retval none
 */

void flag_allocate_sampling(double **rs, double **thetas, double **phis, double **laguweights, double R, int L, int N, enum ssht_methods method){
	
	ssht_allocate_sampling(thetas, phis, L, method);
	flag_spherlaguerre_sampling_allocate(rs, laguweights, N);

}

void flag_deallocate_sampling(double *rs, double *thetas, double *phis, double *laguweights){
	
	free(thetas);
	free(phis);
	flag_spherlaguerre_sampling_deallocate(rs, laguweights);

}

void ssht_allocate_sampling(double **thetas, double **phis, int L, enum ssht_methods method){

	int nphi, ntheta;

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

	  *thetas = (double*)calloc(ntheta, sizeof(double));
	  *phis = (double*)calloc(nphi, sizeof(double));

}

void ssht_sampling(double *thetas, double *phis, int L, enum ssht_methods method){
	
	int t, p, nphi, ntheta;

  	switch (method)
  	{
  	    case MW:
  		    ntheta = ssht_sampling_mw_ntheta(L);
    		nphi = ssht_sampling_mw_nphi(L);
    		for (t=0; t<ntheta; t++)
	      		thetas[t] = ssht_sampling_mw_t2theta(t, L);
	    	for (p=0; p<nphi; p++)
	      		phis[p] = ssht_sampling_mw_p2phi(p, L);
	  		break;
	  	case MWSS:
	  		ntheta = ssht_sampling_mw_ss_ntheta(L);
	    	nphi = ssht_sampling_mw_ss_nphi(L);
	    	for (t=0; t<ntheta; t++)
	      		thetas[t] = ssht_sampling_mw_ss_t2theta(t, L);
		   	for (p=0; p<nphi; p++)
	      		phis[p] = ssht_sampling_mw_ss_p2phi(p, L);
	      	break;
	    case GL:
	    	ntheta = ssht_sampling_gl_ntheta(L);
    		nphi = ssht_sampling_gl_nphi(L);
    		double *weights_unused = (double*)calloc(L,sizeof(double));
    		ssht_sampling_gl_thetas_weights(thetas, weights_unused, L);
    		free(weights_unused);
    		for (p=0; p<nphi; p++)
      			phis[p] = ssht_sampling_gl_p2phi(p, L);
      		break;
      	    case DH:
      		ntheta = ssht_sampling_dh_ntheta(L);
    		nphi = ssht_sampling_dh_nphi(L);
    		for (t=0; t<ntheta; t++)
      			thetas[t] = ssht_sampling_dh_t2theta(t, L);
    		for (p=0; p<nphi; p++)
      			phis[p] = ssht_sampling_dh_p2phi(p, L);
      		break;
	  }
}

void flag_sampling(double *rs, double *thetas, double *phis, double *laguweights, double R, int L, int N, enum ssht_methods method) {
	
	flag_spherlaguerre_sampling(rs, laguweights, R, N);
	ssht_sampling(thetas, phis, L, method);

}
