// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#ifndef FLAG_CORE
#define FLAG_CORE

/*!
 * Allocate FLAG coefficients.
 *
 * \param[out]  flmn Fourier-Laguerre coefficients.
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void flag_allocate_flmn(complex double **flmn, int L, int N);

/*!
 * Allocate sampled field (MW sampling, real).
 *
 * \param[out]  f Sampled dataset (real).
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void flag_allocate_f_real(double **f, int L, int N);

/*!
 * Allocate sampled field (MW sampling, complex).
 *
 * \param[out]  f Sampled dataset (complex).
 * \param[in]  L Angular harmonic band-limit.
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void flag_allocate_f(complex double **f, int L, int N);

/*!
 * Get size of a single layer in MW sampling.
 *
 * \param[in]  L Angular harmonic band-limit.
 * \retval L*(2*L-1) for MW sampling.
 */
int ssht_fr_size(int L);

/*!
 * Get size of a single layer in harmonic space (L^2).
 *
 * \param[in]  L Angular harmonic band-limit.
 * \retval L^2
 */
int ssht_flm_size(int L);

/*!
 * Get size of the full FLAG decomposition.
 *
 * \param[in]  L Angular harmonic band-limit. 
 * \param[in]  N Radial harmonic band-limit.
 * \retval N*L^2
 */
int flag_flmn_size(int L, int N);

/*!
 * Get size of the full dataset for MW sampling.
 *
 * \param[in]  L Angular harmonic band-limit. 
 * \param[in]  N Radial harmonic band-limit.
 * \retval N*L*(2*L-1)
 */
int flag_f_size_mw(int L, int N);

/*!
 * Perform Fourier-Laguerre analysis (MW sampling, complex signal).
 *
 * \param[out] flmn Fourier-Laguerre coefficients.
 * \param[in]  f Input dataset (MW sampling, complex signal)
 * \param[in]  L Angular harmonic band-limit. 
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void flag_analysis(complex double *flmn, const complex double *f, double R, int L, int N);

/*!
 * Perform Fourier-Laguerre synthesis (MW sampling, complex signal).
 *
 * \param[out]  f Input dataset (MW sampling, complex signal)
 * \param[in] flmn Fourier-Laguerre coefficients.
 * \param[in]  L Angular harmonic band-limit. 
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
 void flag_synthesis(complex double *f, const complex double *flmn, const double *nodes, int Nnodes, int L, int N);

/*!
 * Perform Fourier-Laguerre analysis (MW sampling, real signal).
 *
 * \param[out] flmn Fourier-Laguerre coefficients.
 * \param[in]  f Input dataset (MW sampling, real signal)
 * \param[in]  L Angular harmonic band-limit. 
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
void flag_analysis_real(complex double *flmn, const double *f, double R, int L, int N);

/*!
 * Perform Fourier-Laguerre synthesis (MW sampling, real signal).
 *
 * \param[out]  f Input dataset (MW sampling, real signal)
 * \param[in] flmn Fourier-Laguerre coefficients.
 * \param[in]  L Angular harmonic band-limit. 
 * \param[in]  N Radial harmonic band-limit.
 * \retval none
 */
 void flag_synthesis_real(double *f, const complex double *flmn, const double *nodes, int Nnodes, int L, int N);

#endif