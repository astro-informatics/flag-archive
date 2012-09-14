// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flag.h"

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
	for (i=0; i<flag_f_size_mw(L, N); i++){
		f[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
	}
}
	
void flag_random_flmn(complex double *flmn, int L, int N, int seed)
{
	int i;
	srand( time(NULL) );
	for (i=0; i<flag_flmn_size(L, N); i++){
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
	int mapsize = ssht_fr_size_mw(L);
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
	int mapsize = ssht_fr_size_mw(L);
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

void flag_sbesselslag_test_simple(int N, double R, int seed)
{	
	clock_t time_start, time_end;
	int n, i;

	//------- Parameters

	int ell = 0;
	int Nnodes = N;
	int Ninte = 40000;
 	int Nk = N;
 	double kmin = 0.05;
 	double kmax = 0.5;
 	double rmin = 0;
 	double rmax = 5000;


	//------- Input function : f(r) = 1 / r^2

	double *f = (double*)calloc(Nnodes, sizeof(double));
	double *nodes = (double*)calloc(Nnodes, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	flag_spherlaguerre_sampling(nodes, weights, R, N);

	for (n=0; n<Nnodes; n++)
 		f[n] = 1.0 / pow(nodes[n],2.0);

	double *fn = (double*)calloc(N, sizeof(double));
	flag_spherlaguerre_analysis(fn, f, nodes, weights, N);
	double tau = flag_spherlaguerre_tau(R, N);

 	double *kvalues = (double*)calloc(Nk, sizeof(double));
 	for (n=0; n<Nk; n++)
 		kvalues[n] = kmin + (kmax - kmin)*(double)n/(double)Nk;


 	//------- FLK from analytical expression : f_0(k) = sqrt(pi/2) / k

 	double *flk_anal = (double*)calloc(Nk, sizeof(double));
 	for (n=0; n<Nnodes; n++)
 		flk_anal[n] = 1.253314137 / kvalues[n];


 	//------- FLK from discrete integral (naive)

 	double *flk_int = (double*)calloc(Nk, sizeof(double));
 	double r1, r2, f1, f2, h = (rmax - rmin) / Ninte;
 	for (n=0; n<Nnodes; n++){
 		flk_int[n] = 0.0;
		for (i = 0; i < Ninte; i++){
			r1 = rmin + i * h;
		    f1 = sin(r1) / r1;
		    r2 = rmin + (i + 1) * h;
		    f2 = sin(r2) / r2;
	        if(!isnan(f1) && !isinf(f1) && !isnan(f2) && !isinf(f2))
	          flk_int[n] += ((f1 + f2) * h) / (2.0 * kvalues[n]);
		}
	}

	//------- FLK from spherical Laguerre transform

 	double *flk_slag = (double*)calloc(Nk, sizeof(double));
	flag_spherlaguerre2spherbessel(flk_slag, fn, kvalues, Nk, N, ell, tau);


	//-------

	for (n = 0; n < Nk; n++){
		printf("> k = %f : anal = %f / inte = %f / slag = %f\n", kvalues[n], flk_anal[n], flk_int[n], flk_slag[n]);
	}
	for (n = 0; n < Nk; n++){
		printf("> k = %f : slag / anal = %f \n", kvalues[n], flk_slag[n] / flk_anal[n]);
	}
	
	free(f);
	free(fn);
	free(weights);
	free(flk_slag);
	free(flk_anal);
	free(flk_int);
}


void flag_sbesselslag_test_rand(int N, double R, int seed)
{	
	clock_t time_start, time_end;
	double *fn = (double*)calloc(N, sizeof(double));
	int n;
	srand ( time(NULL) );
	for (n=0; n<N; n++)
		fn[n] = rand()/795079784.0;
 	
 	double tau = flag_spherlaguerre_tau(R, N);

 	//-------

 	int Nk = N;
 	double kmin = 0.05;
 	double kmax = 0.5;
 	double Rmin = 0;
 	double Rmax = 2 * R;
 	int ell = 0;
 	int Nnodes = 100;

 	//-------

 	double *flk = (double*)calloc(Nk, sizeof(double));

 	double *kvalues = (double*)calloc(Nk, sizeof(double));
 	for (n=0; n<Nk; n++)
 		kvalues[n] = kmin + (kmax - kmin)*(double)n/(double)Nk;

	time_start = clock();
	//flag_spherlaguerre2spherbessel(flk, fn, kvalues, Nk, N, ell, tau);
	time_end = clock();

	//-------

	double *f = (double*)calloc(Nnodes, sizeof(double));
	double *nodes = (double*)calloc(Nnodes, sizeof(double));
 	for (n=0; n<Nnodes; n++)
 		nodes[n] = Rmin + (Rmax-Rmin) * (double)n / (double)Nnodes;
	flag_spherlaguerre_synthesis(f, fn, nodes, Nnodes, N);

	//-------

 	double *flk_approx = (double*)calloc(Nk, sizeof(double));

	time_start = clock();
	time_end = clock();

	for (n = 0; n < Nk; n++){
		//printf("> k = %f : flk_approx = %f  <->  flk = %f \n", kvalues[n], flk_approx[n], flk[n]);
	}
	
	printf("  - Maximum abs error on reconstruction  : %6.5e\n", 
		maxerr(flk_approx, flk, Nk));

	free(f);
	free(fn);
	free(flk);
	free(flk_approx);
}


int main(int argc, char *argv[]) 
{
	const int NREPEAT = 4;
	const int NSCALE = 5;
	const int L = 4;
	const int N = 8;
	const double R = 420.0;
	const int seed = (int)(10000.0*(double)clock()/(double)CLOCKS_PER_SEC);

	printf("==========================================================\n");
	printf("PARAMETERS : ");
	printf("  L = %i   N = %i   R = %4.1f   seed = %i\n", L, N, R, seed);
	printf("----------------------------------------------------------\n");

	printf("> Testing 1D spherical Bessel transform...\n");
	flag_sbesselslag_test_simple(N, R, seed);
	fflush(NULL);

	printf("==========================================================\n");

	return 0;		
}
