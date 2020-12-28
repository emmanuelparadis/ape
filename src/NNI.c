/* NNI.c    2007-09-05 */

/* Copyright 2007 Vincent Lefort */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "me.h"

//boolean leaf(node *v);
/*edge *siblingEdge(edge *e);
edge *depthFirstTraverse(tree *T, edge *e);
edge *findBottomLeft(edge *e);
edge *topFirstTraverse(tree *T, edge *e);
edge *moveUpRight(edge *e);
double wf(double lambda, double D_LR, double D_LU, double D_LD,
	  double D_RU, double D_RD, double D_DU);*/
/*NNI functions for unweighted OLS topological switches*/

/*fillTableUp fills all the entries in D associated with
  e->head,f->head and those edges g->head above e->head*/
void fillTableUp(edge *e, edge *f, double **A, double **D, tree *T)
{
  edge *g,*h;
  if (T->root == f->tail)
    {
      if (leaf(e->head))
	A[e->head->index][f->head->index] =
	  A[f->head->index][e->head->index] =
	  D[e->head->index2][f->tail->index2];
      else
	{
	  g = e->head->leftEdge;
	  h = e->head->rightEdge;
	  A[e->head->index][f->head->index] =
	    A[f->head->index][e->head->index] =
	    (g->bottomsize*A[f->head->index][g->head->index]
	     + h->bottomsize*A[f->head->index][h->head->index])
	    /e->bottomsize;
	}
    }
  else
    {
      g = f->tail->parentEdge;
      fillTableUp(e,g,A,D,T); /*recursive call*/
      h = siblingEdge(f);
      A[e->head->index][f->head->index] =
	A[f->head->index][e->head->index] =
	(g->topsize*A[e->head->index][g->head->index]
	 + h->bottomsize*A[e->head->index][h->head->index])/f->topsize;
    }
}


void makeOLSAveragesTable(tree *T, double **D, double **A);

double **buildAveragesTable(tree *T, double **D)
{
  int i,j;
  double **A;
  A = (double **) malloc(T->size*sizeof(double *));
  for(i = 0; i < T->size;i++)
    {
      A[i] = (double *) malloc(T->size*sizeof(double));
      for(j=0;j<T->size;j++)
	A[i][j] = 0.0;
    }
  makeOLSAveragesTable(T,D,A);
  return(A);
}

double wf2(double lambda, double D_AD, double D_BC, double D_AC, double D_BD,
	   double D_AB, double D_CD)
{
  double weight;
  weight = 0.5*(lambda*(D_AC + D_BD) + (1 - lambda)*(D_AD + D_BC)
		+ (D_AB + D_CD));
  return(weight);
}

int NNIEdgeTest(edge *e, tree *T, double **A, double *weight)
{
  int a,b,c,d;
  edge *f;
  double *lambda;
  double D_LR, D_LU, D_LD, D_RD, D_RU, D_DU;
  double w1,w2,w0;

  if ((leaf(e->tail)) || (leaf(e->head)))
    return(NONE);
  lambda = (double *)malloc(3*sizeof(double));
  a = e->tail->parentEdge->topsize;
  f = siblingEdge(e);
  b = f->bottomsize;
  c = e->head->leftEdge->bottomsize;
  d = e->head->rightEdge->bottomsize;

  lambda[0] = ((double) b*c + a*d)/((a + b)*(c+d));
  lambda[1] = ((double) b*c + a*d)/((a + c)*(b+d));
  lambda[2] = ((double) c*d + a*b)/((a + d)*(b+c));

  D_LR = A[e->head->leftEdge->head->index][e->head->rightEdge->head->index];
  D_LU = A[e->head->leftEdge->head->index][e->tail->index];
  D_LD = A[e->head->leftEdge->head->index][f->head->index];
  D_RU = A[e->head->rightEdge->head->index][e->tail->index];
  D_RD = A[e->head->rightEdge->head->index][f->head->index];
  D_DU = A[e->tail->index][f->head->index];

  w0 = wf2(lambda[0],D_RU,D_LD,D_LU,D_RD,D_DU,D_LR);
  w1 = wf2(lambda[1],D_RU,D_LD,D_DU,D_LR,D_LU,D_RD);
  w2 = wf2(lambda[2],D_DU,D_LR,D_LU,D_RD,D_RU,D_LD);
  free(lambda);
  if (w0 <= w1)
    {
      if (w0 <= w2) /*w0 <= w1,w2*/
	{
	  *weight = 0.0;
	  return(NONE);
	}
      else /*w2 < w0 <= w1 */
	{
	  *weight = w2 - w0;
/*	  if (verbose)
	    {
	      printf("Swap across %s. ",e->label);
	      printf("Weight dropping by %lf.\n",w0 - w2);
	      printf("New weight should be %lf.\n",T->weight + w2 - w0);
	    }*/
	  return(RIGHT);
	}
    }
  else if (w2 <= w1) /*w2 <= w1 < w0*/
    {
      *weight = w2 - w0;
/*      if (verbose)
	{
	  printf("Swap across %s. ",e->label);
	  printf("Weight dropping by %lf.\n",w0 - w2);
	  printf("New weight should be %lf.\n",T->weight + w2 - w0);
	}*/
      return(RIGHT);
    }
  else /*w1 < w2, w0*/
    {
      *weight = w1 - w0;
/*      if (verbose)
	{
	  printf("Swap across %s. ",e->label);
	  printf("Weight dropping by %lf.\n",w0 - w1);
	  printf("New weight should be %lf.\n",T->weight + w1 - w0);
	}*/
      return(LEFT);
    }
}

int *initPerm(int size);

void NNIupdateAverages(double **A, edge *e, edge *par, edge *skew,
		       edge *swap, edge *fixed, tree *T)
{
  node *v;
  edge *elooper;
  v = e->head;
  /*first, v*/
  A[e->head->index][e->head->index] =
    (swap->bottomsize*
     ((skew->bottomsize*A[skew->head->index][swap->head->index]
       + fixed->bottomsize*A[fixed->head->index][swap->head->index])
      / e->bottomsize) +
     par->topsize*
     ((skew->bottomsize*A[skew->head->index][par->head->index]
       + fixed->bottomsize*A[fixed->head->index][par->head->index])
      / e->bottomsize)
     ) / e->topsize;

  elooper = findBottomLeft(e); /*next, we loop over all the edges
				 which are below e*/
  while (e != elooper)
    {
      A[e->head->index][elooper->head->index] =
	A[elooper->head->index][v->index]
	= (swap->bottomsize*A[elooper->head->index][swap->head->index] +
	   par->topsize*A[elooper->head->index][par->head->index])
	/ e->topsize;
      elooper = depthFirstTraverse(T,elooper);
    }
  elooper = findBottomLeft(swap); /*next we loop over all the edges below and
				    including swap*/
  while (swap != elooper)
  {
    A[e->head->index][elooper->head->index] =
      A[elooper->head->index][e->head->index]
      = (skew->bottomsize * A[elooper->head->index][skew->head->index] +
	 fixed->bottomsize*A[elooper->head->index][fixed->head->index])
      / e->bottomsize;
    elooper = depthFirstTraverse(T,elooper);
  }
  /*now elooper = skew */
  A[e->head->index][elooper->head->index] =
    A[elooper->head->index][e->head->index]
    = (skew->bottomsize * A[elooper->head->index][skew->head->index] +
       fixed->bottomsize* A[elooper->head->index][fixed->head->index])
    / e->bottomsize;

  /*finally, we loop over all the edges in the tree
    on the far side of parEdge*/
  elooper = T->root->leftEdge;
  while ((elooper != swap) && (elooper != e)) /*start a top-first traversal*/
    {
      A[e->head->index][elooper->head->index] =
	A[elooper->head->index][e->head->index]
	= (skew->bottomsize * A[elooper->head->index][skew->head->index]
	   + fixed->bottomsize* A[elooper->head->index][fixed->head->index])
	/ e->bottomsize;
      elooper = topFirstTraverse(T,elooper);
    }

  /*At this point, elooper = par.
    We finish the top-first traversal, excluding the subtree below par*/
  elooper = moveUpRight(par);
  while (NULL != elooper)
    {
      A[e->head->index][elooper->head->index]
	= A[elooper->head->index][e->head->index]
	= (skew->bottomsize * A[elooper->head->index][skew->head->index] +
	   fixed->bottomsize* A[elooper->head->index][fixed->head->index])
	/ e->bottomsize;
      elooper = topFirstTraverse(T,elooper);
    }
}


void NNItopSwitch(tree *T, edge *e, int direction, double **A)
{
  edge *par, *fixed;
  edge *skew, *swap;

/*  if (verbose)
    printf("Branch swap across edge %s.\n",e->label);*/

  if (LEFT == direction)
    swap = e->head->leftEdge;
  else
    swap = e->head->rightEdge;
  skew = siblingEdge(e);
  fixed = siblingEdge(swap);
  par = e->tail->parentEdge;

/*  if (verbose)
    {
      printf("Branch swap: switching edges %s and %s.\n",skew->label,swap->label);
    }*/
  /*perform topological switch by changing f from (u,b) to (v,b)
    and g from (v,c) to (u,c), necessitatates also changing parent fields*/

  swap->tail = e->tail;
  skew->tail = e->head;

  if (LEFT == direction)
    e->head->leftEdge = skew;
  else
    e->head->rightEdge = skew;
  if (skew == e->tail->rightEdge)
    e->tail->rightEdge = swap;
  else
    e->tail->leftEdge = swap;

  /*both topsize and bottomsize change for e, but nowhere else*/

  e->topsize = par->topsize + swap->bottomsize;
  e->bottomsize = fixed->bottomsize + skew->bottomsize;
  NNIupdateAverages(A, e, par, skew, swap, fixed,T);

} /*end NNItopSwitch*/

void reHeapElement(int *p, int *q, double *v, int length, int i);
void pushHeap(int *p, int *q, double *v, int length, int i);
void popHeap(int *p, int *q, double *v, int length, int i);


void NNIRetestEdge(int *p, int *q, edge *e,tree *T, double **avgDistArray,
		double *weights, int *location, int *possibleSwaps)
{
  int tloc;
  tloc = location[e->head->index+1];
  location[e->head->index+1] =
    NNIEdgeTest(e,T,avgDistArray,weights + e->head->index+1);
  if (NONE == location[e->head->index+1])
    {
      if (NONE != tloc)
	popHeap(p,q,weights,(*possibleSwaps)--,q[e->head->index+1]);
    }
  else
    {
      if (NONE == tloc)
	pushHeap(p,q,weights,(*possibleSwaps)++,q[e->head->index+1]);
      else
	reHeapElement(p,q,weights,*possibleSwaps,q[e->head->index+1]);
    }
}

void permInverse(int *p, int *q, int length);

int makeThreshHeap(int *p, int *q, double *v, int arraySize, double thresh);


//void NNI(tree *T, double **avgDistArray, int *count)
void NNI(tree *T, double **avgDistArray, int *count, double **D, int numSpecies)
{
  edge *e, *centerEdge;
  edge **edgeArray;
  int *location;
  int *p,*q;
  int i,j;
  int possibleSwaps;
  double *weights;
  p = initPerm(T->size+1);
  q = initPerm(T->size+1);
  edgeArray = (edge **) malloc((T->size+1)*sizeof(double));
  weights = (double *) malloc((T->size+1)*sizeof(double));
  location = (int *) malloc((T->size+1)*sizeof(int));

  double epsilon = 0.0;
  for (i=0; i<numSpecies; i++)
    for (j=0; j<numSpecies; j++)
      epsilon += D[i][j];
  epsilon = (epsilon / (numSpecies * numSpecies)) * EPSILON;

  for (i=0;i<T->size+1;i++)
    {
      weights[i] = 0.0;
      location[i] = NONE;
    }
  e = findBottomLeft(T->root->leftEdge);
  /* *count = 0; */
  while (NULL != e)
    {
      edgeArray[e->head->index+1] = e;
      location[e->head->index+1] =
	NNIEdgeTest(e,T,avgDistArray,weights + e->head->index + 1);
      e = depthFirstTraverse(T,e);
    }
  possibleSwaps = makeThreshHeap(p,q,weights,T->size+1,0.0);
  permInverse(p,q,T->size+1);
  /*we put the negative values of weights into a heap, indexed by p
    with the minimum value pointed to by p[1]*/
  /*p[i] is index (in edgeArray) of edge with i-th position
    in the heap, q[j] is the position of edge j in the heap */
  while (weights[p[1]] + epsilon < 0)
    {
      centerEdge = edgeArray[p[1]];
      (*count)++;
      T->weight = T->weight + weights[p[1]];
      NNItopSwitch(T,edgeArray[p[1]],location[p[1]],avgDistArray);
      location[p[1]] = NONE;
      weights[p[1]] = 0.0;  /*after the NNI, this edge is in optimal
			      configuration*/
      popHeap(p,q,weights,possibleSwaps--,1);
      /*but we must retest the other four edges*/
      e = centerEdge->head->leftEdge;
      NNIRetestEdge(p,q,e,T,avgDistArray,weights,location,&possibleSwaps);
      e = centerEdge->head->rightEdge;
      NNIRetestEdge(p,q,e,T,avgDistArray,weights,location,&possibleSwaps);
      e = siblingEdge(centerEdge);
      NNIRetestEdge(p,q,e,T,avgDistArray,weights,location,&possibleSwaps);
      e = centerEdge->tail->parentEdge;
      NNIRetestEdge(p,q,e,T,avgDistArray,weights,location,&possibleSwaps);
    }
  free(p);
  free(q);
  free(location);
  free(edgeArray);
}
/*
void NNIwithoutMatrix(tree *T, double **D, int *count)
{
  double **avgDistArray;
  avgDistArray = buildAveragesTable(T,D);
  NNI(T,avgDistArray,count);
}

void NNIWithPartialMatrix(tree *T,double **D,double **A,int *count)
{
  makeOLSAveragesTable(T,D,A);
  NNI(T,A,count);
}
*/
