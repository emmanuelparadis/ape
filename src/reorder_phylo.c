/* reorder_phylo.c       2012-09-03 */

/* Copyright 2008-2012 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>

static int iii;

void foo_reorder(int node, int n, int m, int *e1, int *e2, int *neworder, int *L, int *pos)
{
	int i = node - n - 1, j, k;

/* 'i' is the C index corresponding to 'node' */

	for (j = 0; j < pos[i]; j++) {
		k = L[i + m * j];
		neworder[iii++] = k + 1;
		if (e2[k] > n) /* is it an internal edge? */
			foo_reorder(e2[k], n, m, e1, e2, neworder, L, pos);
	}
}

void bar_reorder(int node, int n, int m, int *e1, int *e2, int *neworder, int *L, int *pos)
{
	int i = node - n - 1, j, k;

	for (j = pos[i] - 1; j >= 0; j--)
		neworder[iii--] = L[i + m * j] + 1;

	for (j = 0; j < pos[i]; j++) {
		k = e2[L[i + m * j]];
		if (k > n)
			bar_reorder(k, n, m, e1, e2, neworder, L, pos);
	}
}

void neworder_phylo(int *n, int *e1, int *e2, int *N, int *neworder, int *order)
/* n: nb of tips
   m: nb of nodes
   N: nb of edges */
{
	int i, j, k, *L, *pos, m = *N - *n + 1, degrmax = *n - m + 1;

/* degrmax is the largest value that a node degree can be */

/* L is a 1-d array storing, for each node, the C indices of the rows of
   the edge matrix where the node is ancestor (i.e., present in the 1st
   column). It is used in the same way than a matrix (which is actually
   a vector) is used in R as a 2-d structure. */

	L = (int*)R_alloc(m * degrmax, sizeof(int));

/* pos gives the position for each 'row' of L, that is the number of elements
   which have already been stored for that 'row'. */

	pos = (int*)R_alloc(m, sizeof(int));
	memset(pos, 0, m * sizeof(int));

/* we now go down along the edge matrix */

	for (i = 0; i < *N; i++) {
		k = e1[i] - *n - 1; /* k is the 'row' index in L corresponding to node e1[i] */
		j = pos[k]; /* the current 'column' position corresponding to k */
		pos[k]++; /* increment in case the same node is found in another row of the edge matrix */
		L[k + m * j] = i;
	}

/* L is now ready: we can start the recursive calls. */
/* We start with the root 'n + 1': its index will be changed into
   the corresponding C index inside the recursive function. */

	switch(*order) {
	case 1 : iii = 0;
		foo_reorder(*n + 1, *n, m, e1, e2, neworder, L, pos);
		break;
	case 2 : iii = *N - 1;
		bar_reorder(*n + 1, *n, m, e1, e2, neworder, L, pos);
		break;
	}
}

#define DO_NODE_PRUNING\
    /* go back down in `edge' to set `neworder' */\
    for (j = 0; j <= i; j++) {\
        /* if find the edge where `node' is */\
        /* the descendant, make as ready */\
        if (edge2[j] == node) ready[j] = 1;\
	if (edge1[j] != node) continue;\
	neworder[nextI] = j + 1;\
	ready[j] = 0; /* mark the edge as done */\
	nextI++;\
    }

void neworder_pruningwise(int *ntip, int *nnode, int *edge1,
			  int *edge2, int *nedge, int *neworder)
{
    int *ready, *Ndegr, i, j, node, nextI, n;

    nextI = *ntip +  *nnode;
    Ndegr = (int*)R_alloc(nextI, sizeof(int));
    memset(Ndegr, 0, nextI*sizeof(int));
    for (i = 0; i < *nedge; i++) (Ndegr[edge1[i] - 1])++;

    ready = (int*)R_alloc(*nedge, sizeof(int));

    /* `ready' indicates whether an edge is ready to be */
    /* collected; only the terminal edges are initially ready */
    for (i = 0; i < *nedge; i++)
	    ready[i] = (edge2[i] <= *ntip) ? 1 : 0;

    /* `n' counts the number of times a node has been seen. */
    /* This algo will work if the tree is in cladewise order, */
    /* so that the nodes of "cherries" will be contiguous in `edge'. */
    n = 0;
    nextI = 0;
    while (nextI < *nedge - Ndegr[*ntip]) {
        for (i = 0; i < *nedge; i++) {
            if (!ready[i]) continue;
	    if (!n) {
	        /* if found an edge ready, initialize `node' and start counting */
	        node = edge1[i];
		n = 1;
	    } else { /* else counting has already started */
	        if (edge1[i] == node) n++;
		else {
		    /* if the node has changed we checked that all edges */
		    /* from `node' have been found */
		    if (n == Ndegr[node - 1]) {
		        DO_NODE_PRUNING
		    }
		    /* in all cases reset `n' and `node' and carry on */
		    node = edge1[i];
		    n = 1;
		}
	    } /* go to the next edge */
	    /* if at the end of `edge', check that we can't do a node */
	    if (n == Ndegr[node - 1]) {
	        DO_NODE_PRUNING
		n = 0;
	    }
        }
    }
    for (i = 0; i < *nedge; i++) {
        if (!ready[i]) continue;
	neworder[nextI] = i + 1;
	nextI++;
    }
}

