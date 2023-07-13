/* dist_nodes.c       2023-07-13 */

/* Copyright 2012-2023 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

/* a and b are integers from the edge matrix, so they start at 1. */
/* NM is a long (64-bit integer) so the result too */
#define DINDEX2(a, b) (a - 1) + NM * (b - 1)

/* The algorithm is pretty simple: the tree must be in cladewise order
   because the edges are visited contiguously. Each edge gives trivially
   one distance, then by moving up along the edge matrix, one finds nodes
   that have already been visited and the distance matrix can be updated. */

SEXP dist_nodes2(SEXP Ntips, SEXP Nnodes, SEXP edge, SEXP edge_length)
{
    /* n: nb of tips, m: nb of nodes, N: nb of edges */

    int i, j, k, a, d, n, m, ROOT, *e1, *e2, N;
    long NM, i1, i2, i3;
    double x, *el, *D;
    SEXP res;

    PROTECT(Ntips = coerceVector(Ntips, INTSXP));
    PROTECT(Nnodes = coerceVector(Nnodes, INTSXP));
    PROTECT(edge = coerceVector(edge, INTSXP));
    PROTECT(edge_length = coerceVector(edge_length, REALSXP));

    n = INTEGER(Ntips)[0];
    m = INTEGER(Nnodes)[0];
    N = LENGTH(edge_length);

    e1 = INTEGER(edge);
    e2 = e1 + N;
    el = REAL(edge_length);

    NM = (long)n + m;

    PROTECT(res = allocMatrix(REALSXP, NM, NM));

    D = REAL(res);
    memset(D, 0, NM * NM * sizeof(double));
    ROOT = e1[0]; d = e2[0]; /* the 2 nodes of the 1st edge */
    i1 = DINDEX2(ROOT, d); i2 = DINDEX2(d, ROOT);
    D[i1] = D[i2] = el[0]; /* the 1st edge gives the 1st distance */

    /* go down along the edge matrix starting at the 2nd edge: */
    for (i = 1; i < N; i++) {
	a = e1[i]; d = e2[i]; x = el[i]; /* get the i-th nodes and branch length */
	i1 = DINDEX2(a, d); i2 = DINDEX2(d, a);
	D[i1] = D[i2] = x;
	/* then go up along the edge matrix from the i-th edge
	   to visit the nodes already visited and update the distances: */
	for (j = i - 1; j >= 0; j--) {
	    k = e2[j];
	    if (k == a) continue;
	    i1 = DINDEX2(k, d); i2 = DINDEX2(d, k); i3 = DINDEX2(a, k);
	    D[i1] = D[i2] = D[i3] + x;
	}
	if (k != ROOT) {
	    i1 = DINDEX2(ROOT, d); i2 = DINDEX2(d, ROOT); i3 = DINDEX2(ROOT, a);
	    D[i1] = D[i2] = D[i3] + x;
	}
    }

    UNPROTECT(5);
    return res;
}
