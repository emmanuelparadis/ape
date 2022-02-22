/* bitsplits.c    2021-12-27 */

/* Copyright 2005-2021 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

/* the following array stores the 8 mask values:
   [0] = 0000 0001
   [1] = 1000 0000
   [2] = 0100 0000
   [3] = 0010 0000
   [4] = 0001 0000
   [5] = 0000 1000
   [6] = 0000 0100
   [7] = 0000 0010
   so that mask81[y % 8] gives the corresponding mask (note that 8 % 8 is 0) */
static const unsigned char mask81[8] = {0x01, 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02};

void OneWiseBitsplits(unsigned char *mat, int nr, int nc, int rest)
{
    /* the following array stores the 8 mask values:
       [0] = 0000 0000
       [1] = 1000 0000
       [2] = 1100 0000
       [3] = 1110 0000
       [4] = 1111 0000
       [5] = 1111 1000
       [6] = 1111 1100
       [7] = 1111 1110
       to set the trailing bits to zero when appropriate */
    const unsigned char trailzeros[8] = {0x00, 0x80, 0xc0, 0xe0, 0xf0, 0xf8, 0xfc, 0xfe};

    int i, j;

    for (i = 0; i < nc; i++) {
	j = nr * i;
	if (mat[j] & mask81[1]) continue;
	while (j < nr * (i + 1)) {
	    mat[j] = ~mat[j];
	    j++;
	}
	if (rest) mat[j - 1] &= trailzeros[rest];
    }
}

static int iii;

void bar_reorder2(int node, int n, int m, int Nedge, int *e, int *neworder, int *L, int *pos)
{
	int i = node - n - 1, j, k;

	for (j = pos[i] - 1; j >= 0; j--)
	    neworder[iii--] = L[i + m * j] + 1;

	for (j = 0; j < pos[i]; j++) {
	    k = e[L[i + m * j] + Nedge];
	    if (k > n)
		bar_reorder2(k, n, m, Nedge, e, neworder, L, pos);
	}
}

#define update_L(x)\
    k = e_reord[i] - Ntip - 1;\
    L[k + Nnode * pos[k]] = x;\
    pos[k]++

SEXP bitsplits_multiPhylo(SEXP x, SEXP n, SEXP nr)
{
    int Ntip, Nnode, Nr, Ntrees, itr, Nc, *e, *e_reord, Nedge, *L, *pos, i, j, k, ispl, *newor, d, inod, y, *rfreq, new_split;
    unsigned char *split, *rmat;
    SEXP mat, freq, ans, EDGE, final_nc;

    PROTECT(x = coerceVector(x, VECSXP));
    PROTECT(n = coerceVector(n, INTSXP)); /* nb of tips */
    PROTECT(nr = coerceVector(nr, INTSXP)); /* nb of rows in the matrix of splits */
    Ntrees = LENGTH(x);
    Ntip = *INTEGER(n);
    Nr = *INTEGER(nr);

    Nc = (Ntip - 3) * Ntrees; /* the maximum number of splits that can be found */

    PROTECT(mat = allocVector(RAWSXP, Nr * Nc));
    PROTECT(freq = allocVector(INTSXP, Nc));
    rmat = RAW(mat);
    rfreq = INTEGER(freq);

    memset(rmat, 0, Nr * Nc * sizeof(unsigned char));

    split = (unsigned char*)R_alloc(Nr, sizeof(unsigned char));

    ispl = 0; /* nb of splits already stored */

    for (itr = 0; itr < Ntrees; itr++) {

	Nnode = *INTEGER(getListElement(VECTOR_ELT(x, itr), "Nnode"));
	if (Nnode == 1) continue;
	PROTECT(EDGE = getListElement(VECTOR_ELT(x, itr), "edge"));
	e = INTEGER(EDGE);
	Nedge = LENGTH(EDGE)/2;

/* L is a 1-d array storing, for each node, the C indices of the rows of
   the edge matrix where the node is ancestor (i.e., present in the 1st
   column). It is used in the same way than a matrix (which is actually
   a vector) is used in R as a 2-d structure. */

	L = (int*)R_alloc(Nnode * Ntip, sizeof(int)); /* safe allocation */

/* pos gives the position for each 'row' of L, that is the number of elements
   which have already been stored for that 'row'. */

	pos = (int*)R_alloc(Nnode, sizeof(int));
	memset(pos, 0, Nnode * sizeof(int));

/* we now go down along the edge matrix */

	for (i = 0; i < Nedge; i++) {
	    k = e[i] - Ntip - 1; /* k is the 'row' index in L corresponding to node e1[i] */
	    j = pos[k]; /* the current 'column' position corresponding to k */
	    pos[k]++; /* increment in case the same node is found in another row of the edge matrix */
	    L[k + Nnode * j] = i;
	}

/* L is now ready: we can start the recursive calls.
   We start with the root 'n + 1': its index will be changed into
   the corresponding C index inside the recursive function. */

	iii = Nedge - 1;
	newor = (int*)R_alloc(Nedge, sizeof(int));
	bar_reorder2(Ntip + 1, Ntip, Nnode, Nedge, e, newor, L, pos);
	e_reord = (int*)R_alloc(2 * Nedge, sizeof(int));
	for (i = 0; i < Nedge; i++) newor[i]--; /* change R indices into C indices */
	for (i = 0; i < Nedge; i++) {
	    e_reord[i] = e[newor[i]];
	    e_reord[i + Nedge] = e[newor[i] + Nedge];
	}
	/* the tree is now reordered */

	/* reallocate L and reinitialize pos */
	L = (int*)R_alloc(Nnode * Ntip, sizeof(int));
	memset(pos, 0, Nnode * sizeof(int));

	for (i = 0; i < Nedge; i++) {
	    memset(split, 0, Nr * sizeof(unsigned char));
	    d = e_reord[i + Nedge];

	    if (d <= Ntip) { /* trivial split from a terminal branch */
		update_L(d);
		continue;
	    }

	    inod = d - Ntip - 1;
	    for (j = 0; j < pos[inod]; j++) {
		y = L[inod + Nnode * j];
		split[(y - 1) / 8] |= mask81[y % 8];
		update_L(y); /* update L */
	    }
	    OneWiseBitsplits(split, Nr, 1, Ntip % 8);
	    new_split = 1;
	    if (itr > 0) {  /* if we are handling the 1st tree, no need to check cause all splits are new */
		j = 0; /* column of rmat */
		k = 0; /* row */
		y = 0; /* number of columns of rmat to shift */
		while (j < ispl) {
		    if (split[k] != rmat[k + y]) { /* the two splits are different so move to the next col of rmat */
			j++;
			k = 0;
			y += Nr;
		    } else k++;

		    if (k == Nr) { /* the two splits are the same, so stop here */
			rfreq[j]++;
			new_split = 0;
			break;
		    }
		}
	    }
	    if (new_split) {
                for (j = 0; j < Nr; j++) rmat[j + ispl * Nr] = split[j];
		rfreq[ispl] = 1;
		ispl++;
	    }
	}
	UNPROTECT(1);
    }
    PROTECT(ans = allocVector(VECSXP, 3));
    PROTECT(final_nc = allocVector(INTSXP, 1));
    INTEGER(final_nc)[0] = ispl;
    SET_VECTOR_ELT(ans, 0, mat);
    SET_VECTOR_ELT(ans, 1, freq);
    SET_VECTOR_ELT(ans, 2, final_nc);
    UNPROTECT(7);
    return ans;
}

int same_splits(unsigned char *x, unsigned char *y, int i, int j, int nr)
{
    int end = i + nr;
    while (i < end) {
	if (x[i] != y[j]) return 0;
	i++;
	j++;
    }
    return 1;
}

SEXP CountBipartitionsFromSplits(SEXP split, SEXP SPLIT)
{
    SEXP FREQ, ans;
    unsigned char *mat, *MAT;
    int i, j, nc, NC, nr, *p, *F;

    PROTECT(split = coerceVector(split, VECSXP));
    PROTECT(SPLIT = coerceVector(SPLIT, VECSXP));
    mat = RAW(getListElement(split, "matsplit"));
    MAT = RAW(getListElement(SPLIT, "matsplit"));

    /* the number of splits in the 1st object: */
    nc = LENGTH(getListElement(split, "freq"));

    /* the split frequencies in the 2nd object: */
    PROTECT(FREQ = getListElement(SPLIT, "freq"));
    F = INTEGER(FREQ);

    /* the number of splits in the 2nd object: */
    NC = LENGTH(FREQ);

    /* the number of rows in the matrix (should be the same in both objects): */
    nr = nrows(getListElement(split, "matsplit"));

    /* create the output */
    PROTECT(ans = allocVector(INTSXP, nc));
    p = INTEGER(ans);
    memset(p, 0, nc * sizeof(int));

    for (i = 0; i < nc; i++) {
	j = 0;
	while (j < NC) {
	    if (same_splits(mat, MAT, nr * i, nr * j, nr)) {
		p[i] = F[j];
		break;
	    }
	    j++;
	}
    }

    UNPROTECT(4);
    return ans;
}
