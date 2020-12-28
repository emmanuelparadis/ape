/* nj.c       2011-10-20 */

/* Copyright 2006-2011 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

double sum_dist_to_i(int n, double *D, int i)
/* returns the sum of all distances D_ij between i and j
   with j = 1...n and j != i */
{
/* we use the fact that the distances are arranged sequentially
   in the lower triangle, e.g. with n = 6 the 15 distances are
   stored as (the C indices are indicated):

           i
     1  2  3  4  5

  2  0
  3  1  5
j 4  2  6  9
  5  3  7 10 12
  6  4  8 11 13 14

  so that we sum the values of the ith column--1st loop--and those of
  (i - 1)th row (labelled 'i')--2nd loop */

	double sum = 0;
	int j, start, end;

	if (i < n) {
		/* the expression below CANNOT be factorized
		   because of the integer operations (it took
		   me a while to find out...) */
		start = n * (i - 1) - i * (i - 1) / 2;
		end = start + n - i;
		for (j = start; j < end; j++) sum += D[j];
	}

	if (i > 1) {
		start = i - 2;
		for (j = 1; j <= i - 1; j++) {
			sum += D[start];
			start += n - j - 1;
		}
	}

	return(sum);
}

SEXP C_nj(SEXP DIST, SEXP N)
{
    double *D, *edge_length, *S, *new_dist, A, B, smallest_S;
    int n, i, j, k, ij, *edge, cur_nod, *otu_label, smallest, OTU1, OTU2, Nedge;
    SEXP phy, E, EL;

    PROTECT(DIST = coerceVector(DIST, REALSXP));
    PROTECT(N = coerceVector(N, INTSXP));
    D = REAL(DIST);
    n = INTEGER(N)[0];

    Nedge = 2 * n - 3;

    PROTECT(phy = allocVector(VECSXP, 2));
    PROTECT(E = allocVector(INTSXP, 2 * Nedge));
    PROTECT(EL = allocVector(REALSXP, Nedge));
    edge = INTEGER(E);
    edge_length = REAL(EL);

    cur_nod = 2 * n - 2;

    S = (double*)R_alloc(n + 1, sizeof(double));
    new_dist = (double*)R_alloc(n * (n - 1) / 2, sizeof(double));
    otu_label = (int*)R_alloc(n + 1, sizeof(int));

    for (i = 1; i <= n; i++) otu_label[i] = i; /* otu_label[0] is not used */

    k = 0;

    while (n > 3) {
	for (i = 1; i <= n; i++) /* S[0] is not used */
	    S[i] = sum_dist_to_i(n, D, i);

	ij = 0;
	smallest_S = 1e50;
	B = n - 2;
	for (i = 1; i < n; i++) {
	    for (j = i + 1; j <= n; j++) {
		A = B * D[ij] - S[i] - S[j];
		if (A < smallest_S) {
		    OTU1 = i;
		    OTU2 = j;
		    smallest_S = A;
		    smallest = ij;
		}
		ij++;
	    }
	}

	edge[k + Nedge] = otu_label[OTU1];
	edge[k + 1 + Nedge] = otu_label[OTU2];
	edge[k] = edge[k + 1] = cur_nod;

	/* get the distances between all OTUs but the 2 selected ones
	   and the latter:
	   a) get the sum for both
	   b) compute the distances for the new OTU */

	A = D[smallest];
	ij = 0;
	for (i = 1; i <= n; i++) {
	    if (i == OTU1 || i == OTU2) continue;
	    new_dist[ij] = (D[give_index(i, OTU1, n)] + /* dist(i, OTU1) */
			    D[give_index(i, OTU2, n)] - /* dist(i, OTU2) */
			    A) / 2;
	    ij++;
	}
	/* compute the branch lengths */
	B = (S[OTU1] - S[OTU2])/B; /* don't need B anymore */
	edge_length[k] = (A + B)/2;
	edge_length[k + 1] = (A - B)/2;

	/* update before the next loop (we are sure that OTU1 < OTU2) */
	if (OTU1 != 1)
	    for (i = OTU1; i > 1; i--) otu_label[i] = otu_label[i - 1];
	if (OTU2 != n)
	    for (i = OTU2; i < n; i++) otu_label[i] = otu_label[i + 1];
	otu_label[1] = cur_nod;

	for (i = 1; i < n; i++) {
	    if (i == OTU1 || i == OTU2) continue;
	    for (j = i + 1; j <= n; j++) {
		if (j == OTU1 || j == OTU2) continue;
		new_dist[ij] = D[DINDEX(i, j)];
		ij++;
	    }
	}

	n--;
	for (i = 0; i < n * (n - 1) / 2; i++) D[i] = new_dist[i];

	cur_nod--;
	k += 2;
    }

    k = 2 * INTEGER(N)[0] - 4; /* 2N - 4 */
    for (i = 0; i < 3; i++) {
	edge[k - i] = cur_nod;
	edge[k - i + Nedge] = otu_label[i + 1];
    }

    edge_length[k] = (D[0] + D[1] - D[2]) / 2;
    k--; /* 2N - 5 */
    edge_length[k] = (D[0] + D[2] - D[1]) / 2;
    k--; /* 2N - 6 */
    edge_length[k] = (D[2] + D[1] - D[0]) / 2;

    SET_VECTOR_ELT(phy, 0, E);
    SET_VECTOR_ELT(phy, 1, EL);

    UNPROTECT(5);
    return phy;
}
