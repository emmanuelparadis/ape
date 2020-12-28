/* bipartition.c    2017-07-28 */

/* Copyright 2005-2017 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

SEXP seq_root2tip(SEXP edge, SEXP nbtip, SEXP nbnode)
{
    int i, j, k, Nedge, *x, *done, dn, sumdone, lt, ROOT, Ntip, Nnode;
    SEXP ans, seqnod, tmp_vec;

    /* The following is needed only if we are not sure
       that the storage mode of `edge' is "integer". */
    PROTECT(edge = coerceVector(edge, INTSXP));
    PROTECT(nbtip = coerceVector(nbtip, INTSXP));
    PROTECT(nbnode = coerceVector(nbnode, INTSXP));
    x = INTEGER(edge); /* copy the pointer */
    Ntip = *INTEGER(nbtip);
    Nnode = *INTEGER(nbnode);
    Nedge = LENGTH(edge)/2;
    ROOT = Ntip + 1;

    PROTECT(ans = allocVector(VECSXP, Ntip));
    PROTECT(seqnod = allocVector(VECSXP, Nnode));

    done = &dn;
    done = (int*)R_alloc(Nnode, sizeof(int));
    for (i = 0; i < Nnode; i++) done[i] = 0;

    tmp_vec = allocVector(INTSXP, 1);
    INTEGER(tmp_vec)[0] = ROOT; /* sure ? */
    SET_VECTOR_ELT(seqnod, 0, tmp_vec);
    sumdone = 0;

    while (sumdone < Nnode) {
        for (i = 0; i < Nnode; i++) { /* loop through all nodes */
	    /* if the vector is not empty and its */
	    /* descendants are not yet found */
	    if (VECTOR_ELT(seqnod, i) == R_NilValue || done[i]) continue;
	    /* look for the descendants in 'edge': */
	    for (j = 0; j < Nedge; j++) {
	        /* skip the terminal edges, we look only for nodes now */
	        if (x[j] - Ntip != i + 1 || x[j + Nedge] <= Ntip) continue;
		/* can now make the sequence from */
		/* the root to the current node */
		lt = LENGTH(VECTOR_ELT(seqnod, i));
		tmp_vec = allocVector(INTSXP, lt + 1);
		for (k = 0; k < lt; k++)
		  INTEGER(tmp_vec)[k] = INTEGER(VECTOR_ELT(seqnod, i))[k];
		INTEGER(tmp_vec)[lt] = x[j + Nedge];
		SET_VECTOR_ELT(seqnod, x[j + Nedge] - Ntip - 1, tmp_vec);
	    }
	    done[i] = 1;
	    sumdone++;
	}
    }

    /* build the sequence from root to tip */
    /* by simply looping through 'edge' */
    for (i = 0; i < Nedge; i++) {
        /* skip the internal edges */
        if (x[i + Nedge] > Ntip) continue;
	lt = LENGTH(VECTOR_ELT(seqnod, x[i] - Ntip - 1));
	tmp_vec = allocVector(INTSXP, lt + 1);
	for (j = 0; j < lt; j++)
	  INTEGER(tmp_vec)[j] = INTEGER(VECTOR_ELT(seqnod, x[i] - Ntip - 1))[j];
	INTEGER(tmp_vec)[lt] = x[i + Nedge];
	SET_VECTOR_ELT(ans, x[i + Nedge] - 1, tmp_vec);
    }

    UNPROTECT(5);
    return ans;
} /* EOF seq_root2tip */
