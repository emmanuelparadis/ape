/* reorder_phylo.c       2021-04-07 */

/* Copyright 2008-2021 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>

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

