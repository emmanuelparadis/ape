/* bitsplits.c    2020-07-19 */

/* Copyright 2005-2020 Emmanuel Paradis */

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

/* #define update_L(x)\ */
/*     k = e[i] - *n - 1;\ */
/*     L[k + *m * pos[k]] = x;\ */
/*     pos[k]++ */

/* void bitsplits_phylo(int *n, int *m, int *e, int *N, int *nr, unsigned char *mat) */
/* /\* n: nb of tips, m: nb of nodes, N: nb of edges, */
/*    nr: number of rows in mat *\/ */
/* { */
/*     int ii, i, j, k, d, y, *L, *pos, inod; */

/*     L = (int*)R_alloc(*n * *m, sizeof(int)); */
/*     pos = (int*)R_alloc(*m, sizeof(int)); */
/*     memset(pos, 0, *m * sizeof(int)); */

/*     ii = 0; */
/*     for (i = 0; i < *N; i++) { */
/* 	d = e[i + *N]; */
/* /\* Rprintf("d = %d\n", d); *\/ */
/* 	if (d <= *n) { /\* trivial split from a terminal branch *\/ */
/* 	    update_L(d); /\* update L *\/ */
/* 	} else { */
/* 	    inod = d - *n - 1; */
/* 	    for (j = 0; j < pos[inod]; j++) { */
/* 		y = L[inod + *m * j]; */
/* /\* Rprintf("\ty = %d\n", y); *\/ */
/* /\* Rprintf("\t\ty / 9 + *nr * ii = %d\n", y / 9 + *nr * ii); *\/ */
/* /\* Rprintf("\t\ty % 8 = %d\n", y % 8); *\/ */
/* 		mat[(y -1) / 8 + *nr * ii] |= mask81[y % 8]; */
/* 		update_L(y); /\* update L *\/ */
/* 	    } */
/* 	    ii++; */
/* 	} */
/*     } */
/*     OneWiseBitsplits(mat, *nr, ii, *n % 8); */
/* } */

/* void CountBipartitionsFromTrees(int *n, int *m, int *e, int *N, int *nr, int *nc, unsigned char *mat, double *freq) */
/* { */
/*     int i, j, k, d, y, *L, *pos, inod; */
/*     unsigned char *split; */

/*     L = (int*)R_alloc(*n * *m, sizeof(int)); */
/*     pos = (int*)R_alloc(*m, sizeof(int)); */
/*     memset(pos, 0, *m * sizeof(int)); */
/*     split = (unsigned char*)R_alloc(*nr, sizeof(unsigned char)); */

/*     for (i = 0; i < *N; i++) { */
/* 	memset(split, 0, *nr * sizeof(unsigned char)); */
/* 	d = e[i + *N]; */
/* 	if (d <= *n) { /\* trivial split from a terminal branch *\/ */
/* 	    update_L(d); */
/* 	} else { */
/* 	    inod = d - *n - 1; */
/* 	    for (j = 0; j < pos[inod]; j++) { */
/* 		y = L[inod + *m * j]; */
/* 		split[(y - 1) / 8] |= mask81[y % 8]; */
/* 		update_L(y); */
/* 	    } */
/* 	} */
/* 	OneWiseBitsplits(split, *nr, 1, *n % 8); */
/* 	j = 0; /\* column of mat *\/ */
/* 	k = 0; /\* row *\/ */
/* 	y = 0; /\* number of columns of mat to shift *\/ */
/* 	while (j < *nc) { */
/* 	    if (split[k] != mat[k + y]) { /\* the two splits are different so move to the next col of mat *\/ */
/* 		j++; */
/* 		k = 0; */
/* 		y += *nr; */
/* 	    } else k++; */
/* 	    if (k == *nr) { /\* the two splits are the same so stop here *\/ */
/* 		freq[j]++; */
/* 		break; */
/* 	    } */
/* 	} */
/*     } */
/* } */

/* #undef update_L */

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

/* Rprintf("Nc = %d\n", Nc); */

    PROTECT(mat = allocVector(RAWSXP, Nr * Nc));
    PROTECT(freq = allocVector(INTSXP, Nc));
    rmat = RAW(mat);
    rfreq = INTEGER(freq);

    memset(rmat, 0, Nr * Nc * sizeof(unsigned char));
    /* memset(rfreq, 0, Nc * sizeof(int)); */

    split = (unsigned char*)R_alloc(Nr, sizeof(unsigned char));

    ispl = 0; /* nb of splits already stored */

    for (itr = 0; itr < Ntrees; itr++) {

/* Rprintf("itr = %d\n", itr); */

	Nnode = *INTEGER(getListElement(VECTOR_ELT(x, itr), "Nnode"));
	PROTECT(EDGE = getListElement(VECTOR_ELT(x, itr), "edge"));
	e = INTEGER(EDGE);
	Nedge = LENGTH(EDGE)/2;

/* Rprintf("Nedge = %d\n", Nedge); */

	/* see explanations in ape/src/reorder_phylo.c */
	L = (int*)R_alloc(Nnode * (Nedge - Ntip + 1), sizeof(int));
	pos = (int*)R_alloc(Nnode, sizeof(int));
	memset(pos, 0, Nnode * sizeof(int));
	for (i = 0; i < Nedge; i++) {
	    k = e[i] - Ntip - 1;
	    j = pos[k];
	    pos[k]++;
	    L[k + Nnode * j] = i;
	}
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
/* Rprintf("itr = %d\ty = %d\n", itr, y); */
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
/* Rprintf("ispl = %d\n", ispl); */
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
