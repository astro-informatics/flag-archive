// FLAG package
// Copyright (C) 2012 
// Boris Leistedt & Jason McEwen

#include "flag.h"

void flag_allocate_flmn(complex double **flmn, int L, int N)
{
	int flmsize = ssht_flm_size(L);
	long totalsize = N*flmsize;
	*flmn = (complex double*)calloc(totalsize, sizeof(complex double));
}

void flag_allocate_f(complex double **f, int L, int N)
{
	int frsize = ssht_fr_size(L);
	long totalsize = N*frsize;
	*f = (complex double*)calloc(totalsize, sizeof(complex double));
}

int ssht_fr_size(int L)
{// In case we want to extend to various sampling schemes
	int mapsize = L*(2*L-1); // MW sampling scheme
	return mapsize;
}

int ssht_flm_size(int L)
{// In case we want to extend to various sampling schemes
	int mapsize = L*L;
	return mapsize;
}

int flag_flmn_size(int L, int N)
{// In case we want to extend to various sampling schemes
	return ssht_flm_size(L)*N;
}

int flag_f_size(int L, int N)
{// In case we want to extend to various sampling schemes
	return ssht_fr_size(L)*N;
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