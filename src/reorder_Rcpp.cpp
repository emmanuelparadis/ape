/* additive.c    2017-07-26 */

/* Copyright 2017 Klaus Schliep */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


static int iii;

void foo_reorderRcpp(int node, int nTips, const IntegerVector & e1,
    const IntegerVector & e2, IntegerVector neworder, const IntegerVector & L,
    const IntegerVector & xi, const IntegerVector & xj)
{
    int i = node - nTips - 1, j, k;
    /* 'i' is the C index corresponding to 'node' */
    for (j = 0; j < xj[i]; j++) {
        k = L[xi[i] + j];
        neworder[iii++] = k + 1;
        if (e2[k] > nTips) /* is it an internal edge? */
        foo_reorderRcpp(e2[k], nTips, e1, e2, neworder, L, xi, xj);
    }
}

void bar_reorderRcpp(int node, int nTips, const IntegerVector & e1,
    const IntegerVector & e2, IntegerVector neworder, const IntegerVector & L,
    const IntegerVector & xi, const IntegerVector & xj)
{
    int i = node - nTips - 1, j, k;

    for (j = xj[i] -1; j >= 0; j--)
        neworder[iii--] = L[xi[i] + j ] + 1;

    for (j = 0; j < xj[i]; j++) {
        k = e2[L[xi[i] + j ]];
        if (k > nTips)
            bar_reorderRcpp(k, nTips, e1, e2, neworder, L, xi, xj);
     }
}


// L is a vector of length number of edges
// not max degree * number of nodes

// [[Rcpp::export]]
IntegerVector reorderRcpp(IntegerMatrix orig, int nTips, int root, int order) {
    IntegerVector e1 = orig( _, 0);
    IntegerVector e2 = orig( _, 1);
    int m = max(e1), k, j;
    int nnode = m - nTips;
//    int root = nTips + 1;
    int n = orig.nrow();
    IntegerVector L(n);
    IntegerVector neworder(n);
    IntegerVector pos(nnode);
    IntegerVector xi(nnode);
    IntegerVector xj(nnode);
    for (int i = 0; i < n; i++) {
        xj[e1[i] - nTips - 1]++;
    }
    for (int i = 1; i < nnode; i++) {
        xi[i] = xi[i-1] + xj[i - 1];
    }
    for (int i = 0; i < n; i++) {
        k = e1[i] - nTips - 1;
        j = pos[k]; /* the current 'column' position corresponding to k */
        L[xi[k] + j] = i;
        pos[k]++;
    }

    switch(order) {
    case 1 : iii = 0;
        foo_reorderRcpp(root, nTips, e1, e2, neworder, L, xi, xj);
        break;
    case 2 : iii = n - 1;
        bar_reorderRcpp(root, nTips, e1, e2, neworder, L, xi, xj);
        break;
    }
    return neworder;
}
