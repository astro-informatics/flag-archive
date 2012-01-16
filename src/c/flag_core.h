
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ssht.h>

#include "flag_spherlaguerre.h"

void flag_analysis(double *flmn, const double *f, int L, int N);

void flag_synthesis(double *f, const double *flmn, int L, int N);
