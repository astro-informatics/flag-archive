
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ssht.h>

#include "flag_sampling.h"
#include "flag_spherlaguerre.h"

/*
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

void ssht_core_forward_real(complex double *flm, double *f, int L, enum ssht_methods method, enum ssht_dl_method_t dl_method){

}
*/
	
void flag_analysis(complex double *flmn, double *f, int L, int N, enum ssht_methods method){

	int n;

	for (n = 0; n < N; n++){
	//	ssht_forward(complex double *flm, complex double *f, int L, int reality, enum ssht_methods method){
	}

}

void flag_synthesis(double *f, complex double *flmn, int L, int N, enum ssht_methods method){
	
}

