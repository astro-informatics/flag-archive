// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flag.h"

/*!
 * Allocate FLAG sampling.
 *
 * \param[out]  rs Radial coordinates.
 * \param[out]  thetas Theta angular coordinates.
 * \param[out]  phis Phi angular coordinates.
 * \param[out]  laguweights Laguerre radial weights for FLAG transform.
 * \param[in]  R Radial boundary / limit.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void flag_allocate_sampling(double **rs, double **thetas, double **phis, double **laguweights, double R, int L, int N)
{
  assert(L > 0);
  assert(N > 1);
  assert(R > 0.0);
	ssht_allocate_sampling(thetas, phis, L);
	flag_allocate_spherlaguerre_sampling(rs, laguweights, N);
}

void ssht_allocate_sampling(double **thetas, double **phis, int L)
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

/*!
 * Compute SSHT MW sampling.
 *
 * \param[out]  thetas Theta angular coordinates.
 * \param[out]  phis Phi angular coordinates.
 * \param[in]  L Angular harmonic band-limit.
 * \retval none
 */
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

/*!
 * Compute FLAG sampling.
 *
 * \param[out]  rs Radial coordinates.
 * \param[out]  thetas Theta angular coordinates.
 * \param[out]  phis Phi angular coordinates.
 * \param[out]  laguweights Laguerre radial weights for FLAG transform.
 * \param[in]  R Radial boundary / limit.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void flag_sampling(double *rs, double *thetas, double *phis, double *laguweights, double R, int L, int N) 
{
  assert(L > 0);
  assert(N > 1);
  assert(R > 0.0);
	flag_spherlaguerre_sampling(rs, laguweights, R, N);
	ssht_sampling_mw(thetas, phis, L);
}
