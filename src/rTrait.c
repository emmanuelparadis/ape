/* rTrait.c       2011-06-25 */

/* Copyright 2010-2011 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>

void C_rTraitCont(int *model, int *Nedge, int *edge1, int *edge2, double *el,
		double *sigma, double *alpha, double *theta, double *x)
{
/* The tree must be in pruningwise order */
	int i;
	double alphaT, M, S;

	switch(*model) {
	case 1 : for (i = *Nedge - 1; i >= 0; i--) {
			GetRNGstate();
			x[edge2[i]] = x[edge1[i]] + sqrt(el[i]) * sigma[i] * norm_rand();
			PutRNGstate();
		}
		break;
	case 2 : for (i = *Nedge - 1; i >= 0; i--) {
			if (alpha[i]) {
				alphaT = alpha[i] * el[i];
				M = exp(-alphaT);
				S = sigma[i] * sqrt((1 - exp(-2 * alphaT))/(2 * alpha[i]));
			} else { /* same than if (alpha[i] == 0) */
				M = 1;
				S = sqrt(el[i]) * sigma[i];
			}
			GetRNGstate();
			x[edge2[i]] = x[edge1[i]] * M +  theta[i] * (1 - M) + S * norm_rand();
			PutRNGstate();
		}
		break;
	}
}
