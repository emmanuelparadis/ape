/* delta_plot.c       2011-06-23 */

/* Copyright 2010-2011 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>

void delta_plot(double *D, int *size, int *nbins, int *counts, double *deltabar)
{
	int x, y, u, v; /* same notation than in Holland et al. 2002 */
	int n = *size, nb = *nbins;
	double dxy, dxu, dxv, dyu, dyv, duv, A, B, C, delta;
	int xy, xu, xv, yu, yv, uv, i;

	for (x = 0; x < n - 3; x++) {
		xy = x*n - x*(x + 1)/2; /* do NOT factorize */
		for (y = x + 1; y < n - 2; y++, xy++) {
			yu = y*n - y*(y + 1)/2; /* do NOT factorize */
			dxy = D[xy];
			xu = xy + 1;
			for (u = y + 1; u < n - 1; u++, xu++, yu++) {
				uv = u*n - u*(u + 1)/2; /* do NOT factorize */
				dxu = D[xu]; dyu = D[yu];
				xv = xu + 1; yv = yu + 1;
				for (v = u + 1; v < n; v++, xv++, yv++, uv++) {
					dxv = D[xv]; dyv = D[yv]; duv = D[uv];
					A = dxv + dyu; B = dxu + dyv; C = dxy + duv;
					if (A == B && B == C) delta = 0; else while (1) {
						if (C <= B && B <= A) {delta = (A - B)/(A - C); break;}
						if (B <= C && C <= A) {delta = (A - C)/(A - B); break;}
						if (A <= C && C <= B) {delta = (B - C)/(B - A); break;}
						if (C <= A && A <= B) {delta = (B - A)/(B - C); break;}
						if (A <= B && B <= C) {delta = (C - B)/(C - A); break;}
						if (B <= A && A <= C) {delta = (C - A)/(C - B); break;}
					}
					/* if (delta == 1) i = nb - 1; else */
					i = delta * nb;
					counts[i] += 1;
					deltabar[x] += delta;
					deltabar[y] += delta;
					deltabar[u] += delta;
					deltabar[v] += delta;
				}
			}
		}
	}
}
