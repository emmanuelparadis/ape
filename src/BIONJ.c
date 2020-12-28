/* BIONJ.c    2012-04-30 */

/* Copyright 2007-2008 Olivier Gascuel, Hoa Sien Cuong,
   R port by Vincent Lefort and Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

/* BIONJ program

   Olivier Gascuel

   GERAD - Montreal- Canada
   olivierg@crt.umontreal.ca

   LIRMM - Montpellier- France
   gascuel@lirmm.fr

   UNIX version, written in C
   by Hoa Sien Cuong (Univ. Montreal) */

#include "me.h"

void Initialize(float **delta, double *X, int n);
void C_bionj(double *X, int *N, int *edge1, int *edge2, double *el);
float Distance(int i, int j, float **delta);
float Variance(int i, int j, float **delta);
int Emptied(int i, float **delta);
float Sum_S(int i, float **delta);
void Compute_sums_Sx(float **delta, int n);
void Best_pair(float **delta, int r, int *a, int *b, int n);
float Agglomerative_criterion(int i, int j, float **delta, int r);
float Branch_length(int a, int b, float **delta, int r);
float Reduction4(int a, float la, int b, float lb, int i, float lamda, float **delta);
float Reduction10(int a, int b, int i, float lamda, float vab, float **delta);
float Lamda(int a, int b, float vab, float **delta, int n, int r);

/* INPUT, OUTPUT, INITIALIZATION

   The lower-half of the delta matrix is occupied by
   dissimilarities. The upper-half of the matrix is
   occupied by variances. The first column
   is initialized as 0; during the algorithm some
   indices are no more used, and the corresponding
   positions in the first column are set to 1. */

/* -- Initialize --

  This function reads an input data and returns the delta matrix

   input: float **delta : delta matrix
          double *X     : distances sent from R as a lower triangle matrix
          int n         : number of taxa

   output: float **delta : delta matrix initialized */

void Initialize(float **delta, double *X, int n)
{
	int i, j; /* matrix line and column indices */
	int k = 0; /* index along X */

	for (i = 1; i < n; i++)
		for (j = i + 1; j <= n; j++)
			delta[i][j] = delta[j][i] = X[k++];

	for (i = 1; i <= n; i++)
		delta[i][i] = delta[i][0] = 0;
}

void C_bionj(double *X, int *N, int *edge1, int *edge2, double *el)
{
	int *a, *b;    /* pair to be agglomerated */
	float **delta; /* delta matrix */
	float la;      /* first taxon branch-length */
	float lb;      /* second taxon branch-length */
	float vab;     /* variance of Dab */
	float lamda = 0.5;

	int r; /* number of subtrees */
	int n; /* number of taxa */
	int i, x, y, curnod, k;

	int *ilab; /* indices of the tips (used as "labels") */

	a = (int*)R_alloc(1, sizeof(int));
	b = (int*)R_alloc(1, sizeof(int));

	n = *N;

	/* Create the delta matrix */
	delta = (float **)R_alloc(n + 1, sizeof(float*));
	for (i = 1; i <= n; i++)
		delta[i] = (float *)R_alloc(n + 1, sizeof(float));

	/* initialise */
	r = n;
	*a = *b = 0;
	Initialize(delta, X, n);

	ilab = (int *)R_alloc(n + 1, sizeof(int));
	for (i = 1; i <= n; i++) ilab[i] = i;

	curnod = 2 * n - 2;
	k = 0;

	while (r > 3) {

		Compute_sums_Sx(delta, n);            /* compute the sum Sx     */
		Best_pair(delta, r, a, b, n);         /* find the best pair by  */
		vab = Variance(*a, *b, delta);        /* minimizing (1)         */
		la = Branch_length(*a, *b, delta, r); /* compute branch-lengths */
		lb = Branch_length(*b, *a, delta, r); /* using formula (2)      */

		lamda = Lamda(*a, *b, vab, delta, n, r); /* compute lambda* using (9)*/

		edge1[k] = edge1[k + 1] = curnod;
		edge2[k] = ilab[*a];
		edge2[k + 1] = ilab[*b];
		el[k] = la;
		el[k + 1] = lb;
		k = k + 2;

		for (i = 1; i <= n; i++) {
			if (Emptied(i,delta) || (i == *a) || (i == *b)) continue;
			if(*a > i) {
				x = *a;
				y = i;
			} else {
				x = i;
				y = *a;
			}

			/* apply reduction formulae 4 and 10 to delta */
			delta[x][y] = Reduction4(*a, la, *b, lb, i, lamda, delta);
			delta[y][x] = Reduction10(*a, *b, i, lamda, vab, delta);
		}

		delta[*b][0] = 1.0; /* make the b line empty */
		ilab[*a] = curnod;
		curnod--;
		r = r - 1;
	}

	/* finalise the three basal edges */
	int last[3];
	i = 1;
	k = 0;
	while (k < 3) {
		if (!Emptied(i, delta)) last[k++] = i;
		i++;
	}

	for (i = 0, k = 2 * n - 4; i < 3; i++, k--) {
		edge1[k] = curnod; /* <- the root at this stage */
		edge2[k] = ilab[last[i]];
	}

	double D[3];
	D[0] = Distance(last[0], last[1], delta);
	D[1] = Distance(last[0], last[2], delta);
	D[2] = Distance(last[1], last[2], delta);

	el[2 * n - 4] = (D[0] + D[1] - D[2])/2;
	el[2 * n - 5] = (D[0] + D[2] - D[1])/2;
	el[2 * n - 6] = (D[2] + D[1] - D[0])/2;
}

/* -- Distance --

   This function retrieves and returns the distance between taxa
   i and j from the delta matrix.

   input: int i         : taxon i
          int j         : taxon j
          float **delta : the delta matrix

   output: float distance : dissimilarity between the two taxa */

float Distance(int i, int j, float **delta)
{
	if (i > j) return(delta[i][j]);
	else return(delta[j][i]);
}

/* -- Variance --

   This function retrieves and returns the variance of the
   distance between i and j, from the delta matrix.

   input: int i         : taxon i
          int j         : taxon j
          float **delta : the delta matrix

   output: float distance : the variance of  Dij */

float Variance(int i, int j, float **delta)
{
	if (i > j) return(delta[j][i]);
	else return(delta[i][j]);
}


/* -- Emptied --

   This function tests if a line is emptied or not.

   input: int i         : subtree (or line) i
          float **delta : the delta matrix

   output: 0 : if not emptied
           1 : if emptied */

int Emptied(int i, float **delta)
{
	return((int)delta[i][0]);
}

/* -- Sum_S --

  This function retrieves the sum Sx from the diagonal
  of the delta matrix

  input: int i         : subtree i
         float **delta : the delta matrix

  output: float delta[i][i] : sum Si */

float Sum_S(int i, float **delta)
{
	return(delta[i][i]);
}

/* -- Compute_sums_Sx --

   This function computes the sums Sx and stores them in the
   diagonal the delta matrix.

   input: float **delta : the delta matrix
          int n         : the number of taxa */


void Compute_sums_Sx(float **delta, int n)
{
	float sum;
	int i, j;

	for (i = 1; i <= n ; i++) {
		if (Emptied(i, delta)) continue;
		sum = 0;
		for (j = 1; j <= n; j++) {
			if (i == j || Emptied(j, delta)) continue;
			sum += Distance(i, j, delta); /* compute the sum Si */
		}
		delta[i][i] = sum;
	}
}


/* -- Best_pair --

   This function finds the best pair to be agglomerated by
   minimizing the agglomerative criterion (1).

   input: float **delta : the delta matrix
          int r         : number of subtrees
          int *a        : contain the first taxon of the pair
          int *b        : contain the second taxon of the pair
          int n         : number of taxa

   output: int *a : the first taxon of the pair
           int *b : the second taxon of the pair */


void Best_pair(float **delta, int r, int *a, int *b, int n)
{
	float Qxy;  /* value of the criterion calculated */
	int x, y;   /* the pair which is tested */
	float Qmin; /* current minimun of the criterion */

	Qmin = 1.0e30;
	for (x = 1; x <= n; x++) {
		if (Emptied(x, delta)) continue;
		for (y = 1; y < x; y++) {
			if (Emptied(y, delta)) continue;
			Qxy = Agglomerative_criterion(x, y, delta, r);
			if (Qxy < Qmin - 0.000001) {
				Qmin = Qxy;
				*a = x;
				*b = y;
			}
		}
	}
}

/* Formulae */

/* Formula (1) */
float Agglomerative_criterion(int i, int j, float **delta, int r)
{
	return((r - 2) * Distance(i, j, delta) -
	       Sum_S(i, delta) - Sum_S(j, delta));
}

/* Formula (2) */
float Branch_length(int a, int b, float **delta, int r)
{
	return(0.5 * (Distance(a, b, delta) +
		      (Sum_S(a, delta) - Sum_S(b, delta))/(r - 2)));
}

/* Formula (4) */
float Reduction4(int a, float la, int b, float lb, int i, float lamda, float **delta)
{
	return(lamda * (Distance(a, i, delta) - la) +
	       (1 - lamda) * (Distance(b, i, delta) - lb));
}

/* Formula (10) */
float Reduction10(int a, int b, int i, float lamda, float vab, float **delta)
{
	return(lamda * Variance(a, i, delta) + (1 - lamda) * Variance(b, i, delta)
	       - lamda * (1 - lamda) * vab);
}

float Lamda(int a, int b, float vab, float **delta, int n, int r)
{
	float lamda = 0.0;
	int i;

	if (vab == 0.0) lamda = 0.5; else {
		for (i = 1; i <= n ; i++) {
			if (a == i || b == i || Emptied(i, delta)) continue;
			lamda += (Variance(b, i, delta) - Variance(a, i, delta));
		}
		lamda = 0.5 + lamda/(2 * (r - 2) * vab); /* Formula (9) */
	}
	if (lamda > 1.0) lamda = 1.0; /*  force 0 < lamda < 1 */
	if (lamda < 0.0) lamda = 0.0;
	return(lamda);
}

