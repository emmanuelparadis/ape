/* matexpo.c       2011-06-23 */

/* Copyright 2007-2011 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <R_ext/Lapack.h>

void mat_expo(double *P, int *nr)
/* This function computes the exponential of a nr x nr matrix */
{
	double *U, *vl, *WR, *Uinv, *WI, *work;
	int i, j, k, l, info, *ipiv, n = *nr, nc = n*n, lw = nc << 1;
	char yes = 'V', no = 'N';

	U = (double *)R_alloc(nc, sizeof(double));
	vl = (double *)R_alloc(n, sizeof(double));
	WR = (double *)R_alloc(n, sizeof(double));
	Uinv = (double *)R_alloc(nc, sizeof(double));
	WI = (double *)R_alloc(n, sizeof(double));
	work = (double *)R_alloc(lw, sizeof(double));

	ipiv = (int *)R_alloc(nc, sizeof(int));

/* The matrix is not symmetric, so we use 'dgeev'.
   We take the real part of the eigenvalues -> WR
   and the right eigenvectors (vr) -> U */
	F77_CALL(dgeev)(&no, &yes, &n, P, &n, WR, WI, vl, &n,
			U, &n, work, &lw, &info);

/* It is not necessary to sort the eigenvalues...
   Copy U -> P */
	memcpy(P, U, nc*sizeof(double));

/* For the inversion, we first make Uinv an identity matrix */
	memset(Uinv, 0, nc*sizeof(double));
	for (i = 0; i < nc; i += n + 1) Uinv[i] = 1;

/* The matrix is not symmetric, so we use 'dgesv'.
   This subroutine puts the result in Uinv (B)
   (P [= U] is erased) */
	F77_CALL(dgesv)(&n, &n, P, &n, ipiv, Uinv, &n, &info);

/* The matrix product of U with the eigenvalues diagonal matrix: */
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			U[j + i*n] *= exp(WR[i]);

/* The second matrix product with U^-1 */
	memset(P, 0, nc*sizeof(double));

	for (k = 0; k < n; k++) {
		for (l = 0; l < n; l++) {
			lw = l + k*n;
			for (i = 0 + n*k, j = l; j < nc; i++, j += n)
				P[lw] += U[j]*Uinv[i];
		}
	}
}
