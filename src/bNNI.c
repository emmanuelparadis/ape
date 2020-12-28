/* bNNI.c    2013-09-26 */

/* Copyright 2007 Vincent Lefort */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "me.h"

/*boolean leaf(node *v);
edge *siblingEdge(edge *e);
edge *depthFirstTraverse(tree *T, edge *e);
edge *findBottomLeft(edge *e);
edge *topFirstTraverse(tree *T, edge *e);
edge *moveUpRight(edge *e);*/

void limitedFillTableUp(edge *e, edge *f, double **A, edge *trigger);
void assignBMEWeights(tree *T, double **A);
//void updateAveragesMatrix(tree *T, double **A, node *v,int direction);
void bNNItopSwitch(tree *T, edge *e, int direction, double **A);
int bNNIEdgeTest(edge *e, tree *T, double **A, double *weight);
void updatePair(double **A, edge *nearEdge, edge *farEdge, node *closer, node *further, double dcoeff, int direction);

int *initPerm(int size);

void reHeapElement(int *p, int *q, double *v, int length, int i);
void pushHeap(int *p, int *q, double *v, int length, int i);
void popHeap(int *p, int *q, double *v, int length, int i);


void bNNIRetestEdge(int *p, int *q, edge *e,tree *T, double **avgDistArray,
		double *weights, int *location, int *possibleSwaps)
{
  int tloc;
  tloc = location[e->head->index+1];
  location[e->head->index+1] =
    bNNIEdgeTest(e,T,avgDistArray,weights + e->head->index+1);
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

int makeThreshHeap(int *p, int *q, double *v, int arraySize, double thresh);

void permInverse(int *p, int *q, int length);

void weighTree(tree *T)
{
  edge *e;
  T->weight = 0;
  for(e = depthFirstTraverse(T,NULL);NULL!=e;e=depthFirstTraverse(T,e))
    T->weight += e->distance;
}

//void bNNI(tree *T, double **avgDistArray, int *count)
void bNNI(tree *T, double **avgDistArray, int *count, double **D, int numSpecies)
{
    edge *e;//, *centerEdge deleted by EP, 2013-09-26, see also below
  edge **edgeArray;
  int *p, *location, *q;
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
/*  if (verbose)
    {
      assignBMEWeights(T,avgDistArray);
      weighTree(T);
    }*/
  e = findBottomLeft(T->root->leftEdge);
  while (NULL != e)
    {
      edgeArray[e->head->index+1] = e;
      location[e->head->index+1] =
	bNNIEdgeTest(e,T,avgDistArray,weights + e->head->index + 1);
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
	/* centerEdge = edgeArray[p[1]]; apparently unused later, deleted by EP, 2013-09-26 */
      (*count)++;
/*      if (verbose)
	{
	  T->weight = T->weight + weights[p[1]];
	  printf("New tree weight is %lf.\n",T->weight);
	}*/
      bNNItopSwitch(T,edgeArray[p[1]],location[p[1]],avgDistArray);
      location[p[1]] = NONE;
      weights[p[1]] = 0.0;  /*after the bNNI, this edge is in optimal
			      configuration*/
      popHeap(p,q,weights,possibleSwaps--,1);
      /*but we must retest the other edges of T*/
      /*CHANGE 2/28/2003 expanding retesting to _all_ edges of T*/
      e = depthFirstTraverse(T,NULL);
      while (NULL != e)
	{
	  bNNIRetestEdge(p,q,e,T,avgDistArray,weights,location,&possibleSwaps);
	  e = depthFirstTraverse(T,e);
	}
    }
  free(p);
  free(q);
  free(location);
  free(edgeArray);
  free(weights);
  assignBMEWeights(T,avgDistArray);
}


/*This function is the meat of the average distance matrix recalculation*/
/*Idea is: we are looking at the subtree rooted at rootEdge.  The subtree
rooted at closer is closer to rootEdge after the NNI, while the subtree
rooted at further is further to rootEdge after the NNI.  direction tells
the direction of the NNI with respect to rootEdge*/
void updateSubTreeAfterNNI(double **A, node *v, edge *rootEdge, node *closer, node *further,
			   double dcoeff, int direction)
{
  edge *sib;
  switch(direction)
    {
    case UP: /*rootEdge is below the center edge of the NNI*/
      /*recursive calls to subtrees, if necessary*/
      if (NULL != rootEdge->head->leftEdge)
	updateSubTreeAfterNNI(A, v, rootEdge->head->leftEdge, closer, further, 0.5*dcoeff,UP);
      if (NULL != rootEdge->head->rightEdge)
	updateSubTreeAfterNNI(A, v, rootEdge->head->rightEdge, closer, further, 0.5*dcoeff,UP);
      updatePair(A, rootEdge, rootEdge, closer, further, dcoeff, UP);
      sib = siblingEdge(v->parentEdge);
      A[rootEdge->head->index][v->index] =
	A[v->index][rootEdge->head->index] =
	0.5*A[rootEdge->head->index][sib->head->index] +
	0.5*A[rootEdge->head->index][v->parentEdge->tail->index];
      break;
    case DOWN: /*rootEdge is above the center edge of the NNI*/
      sib = siblingEdge(rootEdge);
      if (NULL != sib)
	updateSubTreeAfterNNI(A, v, sib, closer, further, 0.5*dcoeff, SKEW);
      if (NULL != rootEdge->tail->parentEdge)
	updateSubTreeAfterNNI(A, v, rootEdge->tail->parentEdge, closer, further, 0.5*dcoeff, DOWN);
      updatePair(A, rootEdge, rootEdge, closer, further, dcoeff, DOWN);
      A[rootEdge->head->index][v->index] =
	A[v->index][rootEdge->head->index] =
	0.5*A[rootEdge->head->index][v->leftEdge->head->index] +
	0.5*A[rootEdge->head->index][v->rightEdge->head->index];
      break;
    case SKEW: /*rootEdge is in subtree skew to v*/
      if (NULL != rootEdge->head->leftEdge)
	updateSubTreeAfterNNI(A, v, rootEdge->head->leftEdge, closer, further, 0.5*dcoeff,SKEW);
      if (NULL != rootEdge->head->rightEdge)
	updateSubTreeAfterNNI(A, v, rootEdge->head->rightEdge, closer, further, 0.5*dcoeff,SKEW);
      updatePair(A, rootEdge, rootEdge, closer, further, dcoeff, UP);
      A[rootEdge->head->index][v->index] =
	A[v->index][rootEdge->head->index] =
	0.5*A[rootEdge->head->index][v->leftEdge->head->index] +
	0.5*A[rootEdge->head->index][v->rightEdge->head->index];
      break;
    }
}

/*swapping across edge whose head is v*/
void bNNIupdateAverages(double **A, node *v, edge *par, edge *skew,
			edge *swap, edge *fixed)
{
  A[v->index][v->index] = 0.25*(A[fixed->head->index][par->head->index] +
				A[fixed->head->index][swap->head->index] +
				A[skew->head->index][par->head->index] +
				A[skew->head->index][swap->head->index]);
  updateSubTreeAfterNNI(A, v, fixed, skew->head, swap->head, 0.25, UP);
  updateSubTreeAfterNNI(A, v, par, swap->head, skew->head, 0.25, DOWN);
  updateSubTreeAfterNNI(A, v, skew, fixed->head, par->head, 0.25, UP);
  updateSubTreeAfterNNI(A, v, swap, par->head, fixed->head, 0.25, SKEW);
}


void bNNItopSwitch(tree *T, edge *e, int direction, double **A)
{
  edge *down, *swap, *fixed;
  node *u, *v;
/*  if (verbose)
    {
      printf("Performing branch swap across edge %s ",e->label);
      printf("with ");
      if (LEFT == direction)
	printf("left ");
      else printf("right ");
      printf("subtree.\n");
    }*/
  down = siblingEdge(e);
  u = e->tail;
  v = e->head;
  if (LEFT == direction)
    {
      swap = e->head->leftEdge;
      fixed = e->head->rightEdge;
      v->leftEdge = down;
    }
  else
    {
      swap = e->head->rightEdge;
      fixed = e->head->leftEdge;
      v->rightEdge = down;
    }
  swap->tail = u;
  down->tail = v;
  if(e->tail->leftEdge == e)
    u->rightEdge = swap;
  else
    u->leftEdge = swap;
  bNNIupdateAverages(A, v, e->tail->parentEdge, down, swap, fixed);
}

double wf5(double D_AD, double D_BC, double D_AC, double D_BD,
	   double D_AB, double D_CD)
{
  double weight;
  weight = 0.25*(D_AC + D_BD + D_AD + D_BC) + 0.5*(D_AB + D_CD);
  return(weight);
}

int bNNIEdgeTest(edge *e, tree *T, double **A, double *weight)
{
  edge *f;
  double D_LR, D_LU, D_LD, D_RD, D_RU, D_DU;
  double w1,w2,w0;
/*  if (verbose)
    printf("Branch swap: testing edge %s.\n",e->label);*/
  if ((leaf(e->tail)) || (leaf(e->head)))
    return(NONE);

  f = siblingEdge(e);

  D_LR = A[e->head->leftEdge->head->index][e->head->rightEdge->head->index];
  D_LU = A[e->head->leftEdge->head->index][e->tail->index];
  D_LD = A[e->head->leftEdge->head->index][f->head->index];
  D_RU = A[e->head->rightEdge->head->index][e->tail->index];
  D_RD = A[e->head->rightEdge->head->index][f->head->index];
  D_DU = A[e->tail->index][f->head->index];

  w0 = wf5(D_RU,D_LD,D_LU,D_RD,D_DU,D_LR); /*weight of current config*/
  w1 = wf5(D_RU,D_LD,D_DU,D_LR,D_LU,D_RD); /*weight with L<->D switch*/
  w2 = wf5(D_DU,D_LR,D_LU,D_RD,D_RU,D_LD); /*weight with R<->D switch*/
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
	      printf("Possible swap across %s. ",e->label);
	      printf("Weight dropping by %lf.\n",w0 - w2);
	      printf("New weight would be %lf.\n",T->weight + w2 - w0);
	    }*/
	  return(RIGHT);
	}
    }
  else if (w2 <= w1) /*w2 <= w1 < w0*/
    {
      *weight = w2 - w0;
/*      if (verbose)
	{
	  printf("Possible swap across %s. ",e->label);
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
	  printf("Possible swap across %s. ",e->label);
	  printf("Weight dropping by %lf.\n",w0 - w1);
	  printf("New weight should be %lf.\n",T->weight + w1 - w0);
	}*/
      return(LEFT);
    }
}

/*limitedFillTableUp fills all the entries in D associated with
  e->head,f->head and those edges g->head above e->head, working
  recursively and stopping when trigger is reached*/
void limitedFillTableUp(edge *e, edge *f, double **A, edge *trigger)
{
  edge *g,*h;
  g = f->tail->parentEdge;
  if (f != trigger)
    limitedFillTableUp(e,g,A,trigger);
  h = siblingEdge(f);
  A[e->head->index][f->head->index] =
    A[f->head->index][e->head->index] =
    0.5*(A[e->head->index][g->head->index] + A[e->head->index][h->head->index]);
}
