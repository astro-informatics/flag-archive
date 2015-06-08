// FLAG package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "flag.h"
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <ssht.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))

double maxerr_cplx(complex double *a, complex double *b, int size)
{
	double value = 0;
	int i;
	for(i = 0; i<size; i++){
		//printf("%6.5e %6.5e %6.5e\n", creal(a[i]-b[i]), cimag(a[i]-b[i]), cabs(a[i]-b[i]));
		value = MAX( cabs( a[i]-b[i] ), value );
	}
	return value;
}

double maxerr(double *a, double *b, int size)
{
	double value = 0;
	int i;
	for(i = 0; i<size; i++){
		value = MAX( abs( a[i]-b[i] ), value );
	}
	return value;
}

double ran2_dp(int idum)
{
  int IM1=2147483563,IM2=2147483399,IMM1=IM1-1,
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
    NTAB=32,NDIV=1+IMM1/NTAB;

  double AM=1./IM1,EPS=1.2e-7,RNMX=1.-EPS;
  int j,k;
  static int iv[32],iy,idum2 = 123456789;
  // N.B. in C static variables are initialised to 0 by default.

  if (idum <= 0) {
    idum= (-idum>1 ? -idum : 1); // max(-idum,1);
    idum2=idum;
    for(j=NTAB+8;j>=1;j--) {
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum=idum+IM1;
      if (j < NTAB) iv[j-1]=idum;
    }
    iy=iv[0];
  }
  k=idum/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum=idum+IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2=idum2+IM2;
  j=1+iy/NDIV;
  iy=iv[j-1]-idum2;
  iv[j-1]=idum;
  if(iy < 1)iy=iy+IMM1;
  return (AM*iy < RNMX ? AM*iy : RNMX); // min(AM*iy,RNMX);

}

void flag_random_f(complex double *f, int L, int N, int seed)
{
	int i;
	srand( time(NULL) );
	for (i=0; i<flag_core_f_size_mw(L, N); i++){
		f[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
	}
}

void flag_random_flmn(complex double *flmn, int L, int N, int seed)
{
	int i;
	srand( time(NULL) );
	for (i=0; i<flag_core_flmn_size(L, N); i++){
		flmn[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);//rand()/795079784.0;
	}
}

void flag_random_flmn_real(complex double *flmn, int L, int N, int seed)
{
	int en, el, m, msign, i, i_op, offset;
	int flmsize = ssht_flm_size(L);

	for (en=0; en<N; en++) {
		offset = en * flmsize;
		for (el=0; el<L; el++) {
    		m = 0;
    		ssht_sampling_elm2ind(&i, el, m);
    		flmn[offset+i] = (2.0*ran2_dp(seed) - 1.0);
    		for (m=1; m<=el; m++) {
      			ssht_sampling_elm2ind(&i, el, m);
      			flmn[offset+i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
      			ssht_sampling_elm2ind(&i_op, el, -m);
      			msign = m & 1;
     			msign = 1 - msign - msign; // (-1)^m
     	 		flmn[offset+i_op] = msign * conj(flmn[offset+i]);
    		}
 		}
 	}
}

void print_f(const complex double *f,int L, int N)
{
	int mapsize = (2*L-1)*L;
	int n, j;
	for(n=0;n<N+1;n++){
		printf("\n -- Layer %i -- \n", n);
		for(j=0;j<mapsize;j++){
			printf(" (%f,%f) ",creal(f[n*mapsize+j]),cimag(f[n*mapsize+j]));
		}
	}
	printf("\n");
}

void print_f_real(const double *f,int L, int N)
{
	int mapsize = (2*L-1)*L;
	int n, j;
	for(n=0;n<N+1;n++){
		printf("\n -- Layer %i -- \n", n);
		for(j=0;j<mapsize;j++){
			printf(" %f ",(f[n*mapsize+j]));
		}
	}
	printf("\n");
}

void print_flmn(const complex double *flmn,int L, int N)
{
	int mapsize = ssht_flm_size(L);
	int n, j;
	for(n=0;n<N;n++){
		printf("\n -- Layer %i -- \n", n);
		for(j=0;j<mapsize;j++){
			printf(" (%f,%f) ",creal(flmn[n*mapsize+j]),cimag(flmn[n*mapsize+j]));
		}
	}
	printf("\n");
}

void flag_sampling_test(int L, int N, double R)
{

	double *rs, *thetas, *phis, *laguweights;

	flag_sampling_allocate(&rs, &thetas, &phis, &laguweights, R, L, N);
	flag_sampling(rs, thetas, phis, laguweights, R, L, N);

	free(rs);
	free(thetas);
	free(phis);
	free(laguweights);

}

void flag_spherlaguerre_quadrature_test(int N)
{
	double *roots = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	const int alpha = ALPHA;

	flag_spherlaguerre_quadrature(roots, weights, N, alpha);
	//int n;
	//for (n=0; n<N+1; n++)
	//	printf("Root %i = %f with weight %f \n",n,roots[n],weights[n]);

}

void flag_spherlaguerre_tau_test(int N)
{
	const double R = 1.0;
	double tau = flag_spherlaguerre_tau(R, N);
	const int alpha = ALPHA;

	double *roots = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	flag_spherlaguerre_quadrature(roots, weights, N, alpha);
	double tau_bis = R / roots[N-1];
	free(roots);
	free(weights);

	//printf("Tau = %4.4e\n",tau);
	//printf("Taubis = %4.4e\n",tau_bis);

	assert(tau != 0);
	assert(cabs(tau_bis-tau) < 1e-14);
}

void flag_spherlaguerre_sampling_test(int N, double R)
{
	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));

	flag_spherlaguerre_sampling(nodes, weights, R, N);
	assert(nodes != NULL);
	assert(weights != NULL);

	//int n;for (n=0; n<N+1; n++)
	//	printf("Node %i = %f \n",n,nodes[n]);
	free(nodes);
	free(weights);
}

void flag_spherlaguerre_cmplx_transform_test(int L, int N, double R, int seed)
{
	clock_t time_start, time_end;
	complex double *f, *fn, *fn_rec;
	int flmsize = ssht_flm_size(L);

	flag_core_allocate_f(&f, L, N);
	flag_core_allocate_flmn(&fn, L, N);
	flag_core_allocate_flmn(&fn_rec, L, N);

	flag_random_flmn(fn, L, N, seed);

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));

 	flag_spherlaguerre_sampling(nodes, weights, R, N);

	time_start = clock();
	flag_spherlaguerre_mapped_synthesis(f, fn, nodes, N, N, flmsize);
	time_end = clock();
	printf("  - Duration of mapped inverse synthesis : %4.4f  seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	flag_spherlaguerre_mapped_analysis(fn_rec, f, nodes, weights, N, flmsize);
	time_end = clock();
	printf("  - Duration of mapped forward analysis  : %4.4f  seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	printf("  - Maximum abs error on reconstruction  : %6.5e\n",
		maxerr_cplx(fn, fn_rec, N));


	free(f);
	free(fn);
	free(fn_rec);
	free(weights);
	free(nodes);
}

void flag_spherlaguerre_transform_test(int N, double R)
{
	clock_t time_start, time_end;
	double *f = (double*)calloc(N, sizeof(double));
	double *fn = (double*)calloc(N, sizeof(double));
	double *fnrec = (double*)calloc(N, sizeof(double));
	int n;

	srand ( time(NULL) );

	for (n=0; n<N; n++){
		fn[n] = rand()/795079784.0;
	}

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));

 	flag_spherlaguerre_sampling(nodes, weights, R, N);

	time_start = clock();
	flag_spherlaguerre_synthesis(f, fn, nodes, N, N);
	time_end = clock();
	printf("  - Duration of 1D inverse synthesis     : %4.0e  seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	flag_spherlaguerre_analysis(fnrec, f, nodes, weights, N);
	time_end = clock();
	printf("  - Duration of 1D forward analysis      : %4.0e  seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	printf("  - Maximum abs error on reconstruction  : %6.5e\n",
		maxerr(fn, fnrec, N));

	/*
	printf("\nTau = %f\n",flag_spherlaguerre_tau(1.0, N));
	for (n=0; n<N+1; n++){
	   printf("%i > node=%f",n, nodes[n]);
	   printf(" weight=%f",weights[n]);
	   printf(" f=%f\n",f[n]);
	}
	for (n=0; n<N; n++){
	   printf(" fn=%f - fnrec=%f",fn[n],fnrec[n]);
	   printf(" d=%2.2e r=%2.2e\n",fn[n]-fnrec[n],fn[n]/fnrec[n]);
	}
	*/


	free(f);
	free(fn);
	free(fnrec);
	free(weights);
	free(nodes);
}

void flag_transform_test(int L, int N, double R, int seed)
{
	complex double *f, *flmn, *flmnrec;
	clock_t time_start, time_end;

	flag_core_allocate_f(&f, L, N);
	flag_core_allocate_flmn(&flmn, L, N);
	flag_core_allocate_flmn(&flmnrec, L, N);

	flag_random_flmn(flmn, L, N, seed);
	//print_flmn(flmn, L, N);

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));

 	flag_spherlaguerre_sampling(nodes, weights, R, N);

    time_start = clock();
	flag_core_synthesis(f, flmn, nodes, N, L, N);
	time_end = clock();
	printf("  - Duration of full 3D synthesis        : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	//print_f(f, L, N);

	time_start = clock();
	flag_core_analysis(flmnrec, f, R, L, N);
	time_end = clock();
	printf("  - Duration of full 3D analysis         : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	//print_flmn(flmnrec, L, N);

	printf("  - Maximum abs error on reconstruction  : %6.5e\n",
		maxerr_cplx(flmn, flmnrec, flag_core_flmn_size(L, N)));

	free(f);
	free(flmn);
	free(flmnrec);
	free(nodes);
	free(weights);
}



void flag_transform_furter_test(int L, int N, double R, int seed)
{
	int spin = 0;
	int verbosity = 0;
	clock_t time_start, time_end, t_for, t_back;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size_mw(L);
	int n, offset_lm, offset_r;
	ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	printf("Sampling...");fflush(NULL);
	flag_spherlaguerre_sampling(nodes, weights, R, N);
	printf("done\n");fflush(NULL);

	complex double *flmn, *flmr, *flmn_rec, *f, *flmr_rec;

	flag_core_allocate_flmn(&flmn, L, N);
	flag_core_allocate_flmn(&flmr, L, N);
	flag_core_allocate_flmn(&flmn_rec, L, N);
	flag_core_allocate_flmn(&flmr_rec, L, N);
	flag_core_allocate_f(&f, L, N);

	flag_random_flmn(flmn, L, N, seed);

	printf("Starting synthesis\n");fflush(NULL);
	time_start = clock();
	flag_spherlaguerre_mapped_synthesis(flmr, flmn, nodes, N, N, flmsize);
	time_end = clock();
	printf("  - Synthesis : SLAG 3D synthes duration : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	time_start = clock();

	flag_spherlaguerre_mapped_analysis(flmn_rec, flmr, nodes, weights, N, flmsize);
	time_end = clock();
	printf("  - Synthesis : SLAG 3D analys duration  : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	printf("  - SLAG maximum absolute error          : %6.5e\n",
		maxerr_cplx(flmn, flmn_rec, flag_core_flmn_size(L, N)));
	//for (n = 0; n < N; n++)
		//for (i = 0; i < flmsize; i++)
			//printf("n=%i - i=%i - (%f,%f) - (%f,%f)\n",n,i,creal(flmn[i+n*flmsize]),cimag(flmn[i+n*flmsize]),creal(flmn_rec[i+n*flmsize]),cimag(flmn_rec[i+n*flmsize]));

	t_for = 0;
	t_back = 0;
	for (n = 0; n < N; n++){
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_inverse_sov_sym(f + offset_r, flmr + offset_lm, L, spin, dl_method, verbosity);
		time_end = clock();
		t_back += (time_end - time_start) / (N);
		time_start = clock();
		ssht_core_mw_forward_sov_conv_sym(flmr_rec + offset_lm, f + offset_r, L, spin, dl_method, verbosity);
		time_end = clock();
		t_for += (time_end - time_start) / (N);
		time_start = clock();
	}
	//for (n = 0; n < N; n++)
		//for (i = 0; i < frsize; i++)
			//printf("n=%i - i=%i - (%f,%f)\n",n,i,creal(f[i+n*frsize]),cimag(f[i+n*frsize]));

	printf("  - Synthesis : SSHT forward av duration : %4.4f seconds\n",
		t_for / (double)CLOCKS_PER_SEC);
	printf("  - Analysis : SSHT backward av duration : %4.4f seconds\n",
		t_back / (double)CLOCKS_PER_SEC);
	printf("  - SSHT maximum absolute error          : %6.5e\n",
		maxerr_cplx(flmr, flmr_rec, flag_core_flmn_size(L, N)));

	free(flmn_rec);
	flag_core_allocate_f(&flmn_rec, L, N);
	free(flmr);
	flag_core_allocate_flmn(&flmr, L, N);

	time_start = clock();
	flag_spherlaguerre_mapped_analysis(flmn_rec, flmr_rec, nodes, weights, N, flmsize);
	time_end = clock();
	printf("  - Analysis : SLAG 3D analysis duration : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	time_start = clock();
	flag_spherlaguerre_mapped_synthesis(flmr, flmn_rec, nodes, N, N, flmsize);
	time_end = clock();
	printf("  - Analysis : SLAG 3D synthesis duration: %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	printf("  - SLAG maximum absolute error          : %6.5e\n",
		maxerr_cplx(flmr, flmr_rec, flag_core_flmn_size(L, N)));
	printf("  - Final FLAG maximum absolute error    : %6.5e\n",
		maxerr_cplx(flmn, flmn_rec, flag_core_flmn_size(L, N)));

	free(flmn);
	free(flmn_rec);
	free(flmr);
	free(flmr_rec);
}

void flag_transform_real_test(int L, int N, double R, int seed)
{
	double *f;
	complex *flmn, *flmnrec;
	clock_t time_start, time_end;

	flag_core_allocate_f_real(&f, L, N);
	flag_core_allocate_flmn(&flmn, L, N);
	flag_core_allocate_flmn(&flmnrec, L, N);

	flag_random_flmn_real(flmn, L, N, seed);
	//print_flmn(flmn, L, N);

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));

 	flag_spherlaguerre_sampling(nodes, weights, R, N);

    time_start = clock();
	flag_core_synthesis_real(f, flmn, nodes, N, L, N);
	time_end = clock();
	printf("  - Duration of full 3D synthesis        : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);

	time_start = clock();
	flag_core_analysis_real(flmnrec, f, R, L, N);
	time_end = clock();
	printf("  - Duration of full 3D analysis         : %4.4f seconds\n",
		(time_end - time_start) / (double)CLOCKS_PER_SEC);
	//print_flmn(flmnrec, L, N);

	printf("  - Maximum abs error on reconstruction  : %6.5e\n",
		maxerr_cplx(flmn, flmnrec, flag_core_flmn_size(L, N)));

	free(f);
	free(flmn);
	free(flmnrec);
	free(nodes);
	free(weights);
}


void flag_transform_performance_test(double R, int NREPEAT, int NSCALE, int seed)
{
	complex double *f, *flmn, *flmnrec;
	clock_t time_start, time_end;
	int sc, repeat;
	double tottime_analysis = 0, tottime_synthesis = 0;
	double accuracy = 0.0;

	int L = 2;
	int N = 2;

	for (sc=0; sc<NSCALE; sc++) {

		L *= 2;
		N *= 2;

		flag_spherlaguerre_tau_test(N);

		flag_core_allocate_flmn(&flmn, L, N);

		double *nodes = (double*)calloc(N, sizeof(double));
		double *weights = (double*)calloc(N, sizeof(double));

		//printf("> Radial sampling...");fflush(NULL);
	 	flag_spherlaguerre_sampling(nodes, weights, R, N);
	 	//printf("done\n");

	 	printf("\n  > R = %4.4f  --  L =  %i N = %i \n", R, L, N);
	 	for (repeat=0; repeat<NREPEAT; repeat++){

	 		//printf("  -> Iteration : %i on %i\n",repeat+1,NREPEAT);

			flag_random_flmn(flmn, L, N, seed);

			flag_core_allocate_f(&f, L, N);

		    time_start = clock();
			flag_core_synthesis(f, flmn, nodes, N, L, N);
			time_end = clock();
			tottime_synthesis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
			//printf("  - Duration for FLAG synthesis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
			flag_core_allocate_flmn(&flmnrec, L, N);

			time_start = clock();
			flag_core_analysis(flmnrec, f, R, L, N);
			time_end = clock();
			tottime_analysis += (time_end - time_start) / (double)CLOCKS_PER_SEC;
			//printf("  - Duration for FLAG analysis   : %4.4f seconds\n", (time_end - time_start) / (double)CLOCKS_PER_SEC);
			accuracy += maxerr_cplx(flmn, flmnrec, flag_core_flmn_size(L, N));
			//printf("  - Max error on reconstruction  : %6.5e\n", maxerr_cplx(flmn, flmnrec, flag_core_flmn_size(L, N)));

			free(f);
			free(flmnrec);

		}

		tottime_synthesis = tottime_synthesis / (double)NREPEAT;
		tottime_analysis = tottime_analysis / (double)NREPEAT;
		accuracy = accuracy / (double)NREPEAT;


		printf("  - Average duration for FLAG synthesis  : %4.4f seconds\n", tottime_synthesis);
		printf("  - Average duration for FLAG analysis   : %4.4f seconds\n", tottime_analysis);
		printf("  - Average max error on reconstruction  : %6.5e\n", accuracy);

		free(flmn);
		free(nodes);
		free(weights);

	}


}


int main(int argc, char *argv[])
{
	//const int NREPEAT = 4;
	//const int NSCALE = 5;
	const int L = 32;
	const int N = 32;
	const double R = 1.0;
	const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);


	double *roots = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	int i;
	//const int alpha = ALPHA;
	//flag_spherlaguerre_quadrature(roots, weights, N, alpha);
    flag_spherlaguerre_sampling(roots, weights, R, N);
	for(i=0;i<N;i++)
		printf("Root[%i] = %f with weight %5.5e\n",i,roots[i],weights[i]);
	free(roots);
	free(weights);


	printf("==========================================================\n");
	printf("PARAMETERS : ");
	printf("  L = %i   N = %i   R = %4.1f   seed = %i\n", L, N, R, seed);
	printf("----------------------------------------------------------\n");
	printf("> Testing Laguerre quadrature...");
	flag_spherlaguerre_quadrature_test(N);
	printf("OK\n");

	printf("> Testing Laguerre tau calculation...");
	flag_spherlaguerre_tau_test(N);
	printf("OK\n");

	printf("> Testing Laguerre sampling scheme...");
	flag_spherlaguerre_sampling_test(N, R);
	printf("OK\n");

	printf("----------------------------------------------------------\n");

	printf("> Testing 1D Laguerre transform...\n");
	flag_spherlaguerre_transform_test(N, R);
	fflush(NULL);

	printf("----------------------------------------------------------\n");

	printf("> Testing cmplx mapped Laguerre transform...\n");
	flag_spherlaguerre_cmplx_transform_test(L, N, R, seed);
	fflush(NULL);

	printf("----------------------------------------------------------\n");

	printf("> Testing FLAG sampling scheme...");
	flag_sampling_test(L, N, R);
	printf("OK\n");

	printf("> Testing FLAG transform...\n");
	flag_transform_test(L, N, R, seed);
	fflush(NULL);

	printf("> Testing REAL FLAG transform...\n");
	flag_transform_real_test(L, N, R, seed);
	fflush(NULL);

	printf("----------------------------------------------------------\n");

	printf("> Testing FLAG in further details...\n");
	flag_transform_furter_test(L, N, R, seed);
	fflush(NULL);

	printf("==========================================================\n");

	//printf("> FLAG transform : performance and accuracy tests\n");
	//flag_transform_performance_test(R, NREPEAT, NSCALE, seed);
	//fflush(NULL);

	//printf("==========================================================\n");

	return 0;
}
