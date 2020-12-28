/* ape.h    2014-01-02 */

/* Copyright 2011-2014 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

#define DINDEX(i, j) n*(i - 1) - i*(i - 1)/2 + j - i - 1

/* in ape.c */
int give_index(int i, int j, int n);
SEXP getListElement(SEXP list, char *str);

/* in njs.c */
void choosePair(double* D, int n, double* R, int* s, int* sw, int* x, int* y, int fS);
double cnxy(int x, int y, int n, double* D);
int mxy(int x,int y, int n, double* D);
double nxy(int x, int y, int n, double* D);
int cxy(int x, int y, int n, double* D);

/* in triangMtd.c */
void C_triangMtd(double* d, int* np, int* ed1, int* ed2, double* edLen);
int * getPathBetween(int x, int y, int n, int* ed1, int* ed2, int numEdges);
int give_indexx(int i, int j, int n); /* a variant of the above */

