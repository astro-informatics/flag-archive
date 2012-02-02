
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <flag.h>

/*void assert_double(double a, double b, String code){
	if( a != b ){
		FLAG_ERROR_GENERIC(code)
	}
}*/

void flag_sampling_test(int L, int N){

	const double R = 1.0;
	double *rs, *thetas, *phis, *laguweights;
	int i, n;

	for( i=MW; i <= DH; i++ ) {

	   flag_allocate_sampling(&rs, &thetas, &phis, &laguweights, R, L, N, MW);
	   flag_sampling(rs, thetas, phis, laguweights, R, L, N, MW);
	   //for (n=0; n<N; n++)
		//printf("r=%f t=%f p=%f \n",rs[n], thetas[n], phis[n]);
	   flag_deallocate_sampling(rs, thetas, phis, laguweights);

	}
	
}

void flag_spherlaguerre_quadrature_test(int N){

	double *roots = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));

	flag_spherlaguerre_quadrature(roots, weights, N);

	/*
	int n;
	for (n=0; n<N; n++)
		printf("Root %i = %f with weight %f \n",n,roots[n],weights[n]);
	*/

}

void flag_spherlaguerre_tau_test(int N){
	
	const double R = 1.0;
	double tau = flag_spherlaguerre_tau(R,N);
	//printf("Tau = %f",tau);

}

void flag_spherlaguerre_sampling_test(int N){

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));
	const double R = 1.0;

	flag_spherlaguerre_sampling(nodes, weights, R, N);

	
	int n;
	//for (n=0; n<N; n++)
	//	printf("Node %i = %f \n",n,nodes[n]);
	

	free(nodes);

}

void flag_spherlaguerre_transform_test(int N){
	
	double *f = (double*)calloc(N, sizeof(double));
	double *fn = (double*)calloc(N, sizeof(double));
	double *frec = (double*)calloc(N, sizeof(double));
	const double R = 1.0;
	int n;

	srand ( time(NULL) );

	for (n=0; n<N; n++){
		f[n] = rand()/795079784.0;
	}

	double *nodes = (double*)calloc(N, sizeof(double));
	double *weights = (double*)calloc(N, sizeof(double));

 	flag_spherlaguerre_sampling(nodes, weights, R, N);
		
	flag_spherlaguerre_analysis(fn, f, nodes, weights, N);

	flag_spherlaguerre_synthesis(frec, fn, nodes, N);
	
	//printf("\nTau = %f\n",flag_spherlaguerre_tau(1.0, N));

	for (n=0; n<N; n++){
	   /*
	   printf("\nnodes[%i] = %f",n,nodes[n]);
	   printf("\nweights[%i] = %f",n,weights[n]);
	   printf("\nfn[%i] = %f",n,fn[n]);
	   */
	   //printf("\nf[%i] = %f - frec[%i] = %f",n,f[n],n,frec[n]);
	}

	free(f);
	free(fn);
	free(frec);
	free(weights);
	free(nodes);

}

int main(int argc, char *argv[]) {

	const int L = 10;
	const int N = 10;
	const double R = 1.0;

	printf("> Testing Laguerre quadrature...");
	flag_spherlaguerre_quadrature_test(N);
	printf("OK\n");

	printf("> Testing Laguerre tau calculation...");
	flag_spherlaguerre_tau_test(N);
	printf("OK\n");

	printf("> Testing Laguerre sampling scheme...");
	flag_spherlaguerre_sampling_test(N);
	printf("OK\n");

	printf("> Testing Laguerre transform...");
	flag_spherlaguerre_transform_test(N);
	printf("OK\n");

	printf("> Testing FLAG sampling schemes...");
	flag_sampling_test(L, N);
	printf("OK\n");
	
	return 0;		
}
