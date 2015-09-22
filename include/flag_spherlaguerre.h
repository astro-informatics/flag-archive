
#ifndef FLAG_SPHERLAGUERRE
#define FLAG_SPHERLAGUERRE

#include <complex.h>

/*!
 * Compute Gauss-Laguerre quadrature (nodes and weights).
 *
 * \param[out]  roots Gauss-Laguerre nodes.
 * \param[out]  weights Gauss-Laguerre weights.
 * \param[in]  N Harmonic band-limit.
 * \retval none
 */
void flag_spherlaguerre_quadrature(double *roots, double *weights, int N, int alpha);

/*!
 * Compute spherical Laguerre scaling factor tau.
 *
 * \param[in]  R Radial limit / boundary.
 * \param[in]  N Harmonic band-limit.
 * \retval tau Scaling factor for the SLAG sampling.
 */
double flag_spherlaguerre_tau(double R, int N);

/*!
 * Compute spherical Laguerre maximum R.
 *
 * \param[in]  N Harmonic band-limit.
 * \retval R Radial limit / boundary.
 */
double flag_spherlaguerre_Rmax(int N);

/*!
 * Compute spherical Laguerre sampling scheme.
 *
 * \param[out]  nodes Nodes of the sampling.
 * \param[out]  weights Weights for the SLAG quadrature.
 * \param[in]  R Radial limit / boundary.
 * \param[in]  N Harmonic band-limit.
 * \retval none
 */
void flag_spherlaguerre_sampling(double *nodes, double *weights, double tau, int N);

void flag_spherbessel_sampling(double *nodes, double *weights, double R, int N);

/*!
 * Perform spherical Laguerre analysis.
 *
 * \param[out]  fn SLAG coefficients.
 * \param[in]  f Input dataset.
 * \param[in]  nodes Nodes of the sampling.
 * \param[in]  weights Weights for the SLAG quadrature.
 * \param[in]  N Harmonic band-limit.
 * \retval none
 */
void flag_spherlaguerre_analysis(double *fn, const double *f, const double *nodes, const double *weights, double tau, int N);

/*!
 * Perform spherical Laguerre synthesis.
 *
 * \param[out]  f Synthesised dataset.
 * \param[in]  fn Input SLAG coefficients.
 * \param[in]  nodes Nodes of the sampling.
 * \param[in]  N Harmonic band-limit.
 * \retval none
 */
void flag_spherlaguerre_synthesis(double *f, const double *fn, const double *nodes, int Nnodes, double tau, int N);

void flag_spherlaguerre_synthesis_gen(double *f, const double *fn, const double *nodes, int Nnodes, double tau, int N, int alpha);


/*!
 * Allocate spherical Laguerre sampling scheme.
 *
 * \param[out]  nodes Nodes of the sampling.
 * \param[out]  weights Weights for the SLAG quadrature.
 * \param[in]  N Harmonic band-limit.
 * \retval none
 */
void flag_spherlaguerre_allocate_sampling(double **nodes, double **weights, int N);

/*!
 * Perform spherical Laguerre analysis.
 * 3D mapped version - suitable for FLAG transform.
 *
 * \param[out]  fn SLAG coefficients.
 * \param[in]  f Input dataset.
 * \param[in]  nodes Nodes of the sampling.
 * \param[in]  weights Weights for the SLAG quadrature.
 * \param[in]  mapsize Size of each layer (L^2 in FLAG).
 * \param[in]  N Harmonic band-limit.
 * \retval none
 */
void flag_spherlaguerre_mapped_analysis(complex double *fn, const complex double *f, const double *weights, const double *nodes, double tau, int N, int mapsize);

/*!
 * Perform spherical Laguerre synthesis.
 * 3D mapped version - suitable for FLAG transform.
 *
 * \param[out]  f Synthesised dataset.
 * \param[in]  fn Input SLAG coefficients.
 * \param[in]  nodes Nodes of the sampling.
 * \param[in]  mapsize Size of each layer (L^2 in FLAG).
 * \param[in]  N Harmonic band-limit.
 * \retval none
 */
void flag_spherlaguerre_mapped_synthesis(complex double *f, const complex double *fn, const double *nodes, int Nnodes, double tau, int N, int mapsize);


/*!
 * Compute the spherical-Laguerre basis functions on a grid of radii
 *
 * \param[out]  KN Synthesised basis function.
 * \param[in]  N SLAG band-limit.
 * \param[in]  nodes Radii where the function must be calculated.
 * \param[in]  Nnodes Number of radii.
 * \retval none
 */
void flag_spherlaguerre_basis(double *KN, const int N, const double *nodes, int Nnodes, double tau);

#endif
