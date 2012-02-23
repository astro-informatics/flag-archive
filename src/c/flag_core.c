// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flag.h"

double ran2_dp(int idum) {

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
	for (i=0; i<N*ssht_fr_size(L); i++){
		f[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
	}
}
	
void flag_random_flmn(complex double *flmn, int L, int N, int seed)
{
	int i;
	srand( time(NULL) );
	for (i=0; i<N*ssht_flm_size(L); i++){
		flmn[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);//rand()/795079784.0;
	}
}

void flag_allocate_flmn(complex double **flmn, int L, int N)
{
	int flmsize = ssht_flm_size(L);
	//printf("\n FLMSIZE = %i \n",flmsize);
	long totalsize = N*flmsize;
	*flmn = (complex double*)calloc(totalsize, sizeof(complex double));
}

void flag_allocate_f(complex double **f, int L, int N)
{
	int frsize = ssht_fr_size(L);
	//printf("\n FRSIZE = %i \n",frsize);
	long totalsize = N*frsize;
	*f = (complex double*)calloc(totalsize, sizeof(complex double));
}

int ssht_fr_size(int L)
{
	// In case we want to extend to various sampling schemes
	int mapsize = L*(2*L-1); // MW sampling scheme
	return mapsize;
}

int ssht_flm_size(int L)
{
	// In case we want to extend to various sampling schemes
	int mapsize = L*L;
	return mapsize;
}

void flag_analysis(complex double *flmn, 
		const complex double *f, 
		int L, int N)
{

	int spin = 0;
	int verbosity = 0;
	int i, n, offset_lm, offset_r;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size(L);
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	complex double *flmr;
	flag_allocate_flmn(&flmr, L, N);

	for (n = 0; n < N; n++){
		//printf("> Analysis: layer %i on %i\n",n+1,N);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_forward_sov_conv_sym(flmr + offset_lm, f + offset_r, L, spin, dl_method, verbosity);
	}


	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	flag_spherlaguerre_quadrature(nodes, weights, N);
	flag_mapped_spherlaguerre_analysis(flmn, flmr, nodes, weights, flmsize, N);
	free(nodes);
	free(weights);
    free(flmr);

}

void flag_synthesis(complex double *f, 
		const complex double *flmn, 
		int L, int N)
{
	
	int spin = 0;
	int verbosity = 0;
	int n, i, offset_lm, offset_r;
	int flmsize = ssht_flm_size(L);
	int frsize = ssht_fr_size(L);
	ssht_dl_method_t dl_method = SSHT_DL_RISBO;

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	flag_spherlaguerre_quadrature(nodes, weights, N);

	complex double *flmr;
	flag_allocate_flmn(&flmr, L, N);
	flag_mapped_spherlaguerre_synthesis(flmr, flmn, nodes, flmsize, N);
	free(nodes);
	free(weights);

	for (n = 0; n < N; n++){
		//printf("> Synthesis: layer %i on %i\n",n+1,N);
		offset_lm = n * flmsize;
		offset_r = n * frsize;
		ssht_core_mw_inverse_sov_sym(f + offset_r, flmr + offset_lm, L, spin, dl_method, verbosity);
	}

    free(flmr);

}

void flag_analysis_real(complex double *flmn, 
		const double *f, 
		int L, int N)
{
			
}

void flag_synthesis_real(double *f, 
		const complex double *flmn, 
		int L, int N)
{

}


/* In case we want to extend to various sampling schemes
void ssht_core_forward(complex double *flm, complex double *f, int L, enum ssht_methods method, enum ssht_dl_method_t dl_method){
	
	int spin = 0;
	int verbosity = 0;

	// TODO : NORTH / SOUTH POLE PROBLEM
	// TODO : SPIN PROBLEM

	switch (method)
  	{
  		case MW:
  			if (reality)
				ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, dl_method, verbosity);
      		else
				ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, dl_method, verbosity);
    		break;
	  	case MWSS:
	  		if (reality)
				ssht_core_mw_forward_sov_conv_sym_ss_real(flm, f, L, dl_method, verbosity);   
      		else
				ssht_core_mw_forward_sov_conv_sym_ss(flm, f, L, spin,  dl_method, verbosity);   
    		break;
	    case GL:
	    	if (reality)
      			ssht_core_gl_forward_sov_real(flm, f, L, verbosity);    
    		else
      			ssht_core_gl_forward_sov(flm, f, L, spin, verbosity);
    		break;
      	case DH:
      		if (reality)
      			ssht_core_dh_forward_sov_real(flm, f, L, verbosity);    
    		else
      			ssht_core_dh_forward_sov(flm, f, L, spin, verbosity);    
    		break;
      	case HPX: //TODO
      		break;
	  	default: //TODO
	  		break;
	  }

}

*/