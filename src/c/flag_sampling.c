// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flag.h"

/*!  
 * Compute the exact Fourier-Laguerre sampling.
 *
 * \param[out] f Sampling on the 3-ball.
 * \param[in] L Harmonic band-limit.
 * \param[in] N Laguerre band-limit.
 * \retval none
 */

void flag_allocate_sampling(double **rs, double **thetas, double **phis, double **laguweights, double R, int L, int N)
{
	ssht_allocate_sampling(thetas, phis, L);
	flag_allocate_spherlaguerre_sampling(rs, laguweights, N);
}

void ssht_allocate_sampling(double **thetas, double **phis, int L)
{

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

void ssht_sampling(double *thetas, double *phis, int L)
{
	
	int t, p, nphi, ntheta;

	ntheta = ssht_sampling_mw_ntheta(L);
  nphi = ssht_sampling_mw_nphi(L);
  for (t=0; t<ntheta; t++)
	 	thetas[t] = ssht_sampling_mw_t2theta(t, L);
	for (p=0; p<nphi; p++)
	  phis[p] = ssht_sampling_mw_p2phi(p, L);

  	/* FOR FUTURE IMPROVEMENTS // multi-scheme support
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
	  */

}

void flag_sampling(double *rs, double *thetas, double *phis, double *laguweights, double R, int L, int N) 
{
	flag_spherlaguerre_sampling(rs, laguweights, R, N);
	ssht_sampling(thetas, phis, L);
}
