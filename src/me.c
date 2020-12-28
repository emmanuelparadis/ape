/* me.c    2019-03-26 */

/* Copyright 2007-2008 Olivier Gascuel, Rick Desper,
   R port by Vincent Lefort and Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "me.h"

//functions from me_balanced.c
tree *BMEaddSpecies(tree *T, node *v, double **D, double **A);
void assignBMEWeights(tree *T, double **A);
void makeBMEAveragesTable(tree *T, double **D, double **A);
//functions from me_ols.c
tree *GMEaddSpecies(tree *T, node *v, double **D, double **A);
void assignOLSWeights(tree *T, double **A);
void makeOLSAveragesTable(tree *T, double **D, double **A);
//functions from bNNI.c
void bNNI(tree *T, double **avgDistArray, int *count, double **D, int numSpecies);
//functions from NNI.c
void NNI(tree *T, double **avgDistArray, int *count, double **D, int numSpecies);
//functions from SPR.c
void SPR(tree *T, double **D, double **A, int *count);
//functions from TBR.c
//void TBR(tree *T, double **D, double **A);


void me_b(double *X, int *N, int *labels,
	  int *nni, int *spr, int *tbr,
	  int *edge1, int *edge2, double *el)
{
  double **D, **A;
  set *species, *slooper;
  node *addNode;
  tree *T;
  int n, nniCount;

  n = *N;
  T = NULL;
  nniCount = 0;
  species = (set *) malloc(sizeof(set));
  species->firstNode = NULL;
  species->secondNode = NULL;
  D = loadMatrix(X, labels, n, species);
  A = initDoubleMatrix(2*n - 2);

  for(slooper = species; NULL != slooper; slooper = slooper->secondNode)
  {
    addNode = copyNode(slooper->firstNode);
    T = BMEaddSpecies(T, addNode, D, A);
  }
  // Compute bNNI
  if (*nni) bNNI(T, A, &nniCount, D, n);
  assignBMEWeights(T,A);

  if (*spr) SPR(T, D, A, &nniCount);
  if (*tbr) Rprintf("argument tbr was ignored: TBR not performed\n"); //TBR(T, D, A);

  tree2phylo(T, edge1, edge2, el, labels, n);

  freeMatrix(D,n);
  freeMatrix(A,2*n - 2);
  freeSet(species);
  freeTree(T);
  T = NULL;
}

void me_o(double *X, int *N, int *labels, int *nni,
	  int *edge1, int *edge2, double *el)
{
  double **D, **A;
  set *species, *slooper;
  node *addNode;
  tree *T;
  int n, nniCount;

  n = *N;
  T = NULL;
  nniCount = 0;
  species = (set *) malloc(sizeof(set));
  species->firstNode = NULL;
  species->secondNode = NULL;

  D = loadMatrix (X, labels, n, species);
  A = initDoubleMatrix(2 * n - 2);

  for(slooper = species; NULL != slooper; slooper = slooper->secondNode)
  {
    addNode = copyNode(slooper->firstNode);
    T = GMEaddSpecies(T,addNode,D,A);
  }
  makeOLSAveragesTable(T,D,A);
  // Compute NNI
  if (*nni)
    NNI(T,A,&nniCount,D,n);
  assignOLSWeights(T,A);

  tree2phylo(T, edge1, edge2, el, labels, n);

  freeMatrix(D,n);
  freeMatrix(A,2*n - 2);
  freeSet(species);
  freeTree(T);
  T = NULL;
}

/*

  -- MATRIX FUNCTIONS --

*/

double **initDoubleMatrix(int d)
{
  int i,j;
  double **A;
  A = (double **) malloc(d*sizeof(double *));
  for(i=0;i<d;i++)
    {
      A[i] = (double *) malloc(d*sizeof(double));
      for(j=0;j<d;j++)
	A[i][j] = 0.0;
    }
  return(A);
}

//double **loadMatrix (double *X, char **labels, int n, set *S)
double **loadMatrix (double *X, int *labels, int n, set *S)
{
//  char nextString[MAX_LABEL_LENGTH];
  node *v;
  double **table;
  int i, j, a, b;

  table = (double **) calloc(n,sizeof(double *));
  for(i=0; i<n; i++)
    table[i] = (double *) calloc(n,sizeof(double));

  for(i=0; i<n; i++)
    {
//      strncpy (nextString, labels[i], MAX_LABEL_LENGTH);
//      ReplaceForbiddenChars (nextString, '_');
//      v = makeNewNode(nextString,-1);
      v = makeNewNode(labels[i], -1);
      v->index2 = i;
      S = addToSet(v,S);
      for (j=i; j<n; j++) {
        a=i+1;
        b=j+1;
        table[j][i] = X[XINDEX(a,b)];
        table[i][j] = X[XINDEX(a,b)];
        if (i==j)
          table[i][j] = 0;
      }
    }
  return (table);
}

/*

  -- GRAPH FUNCTIONS --

*/

set *addToSet(node *v, set *X)
{
  if (NULL == X)
    {
      X = (set *) malloc(sizeof(set));
      X->firstNode = v;
      X->secondNode = NULL;
    }
  else if (NULL == X->firstNode)
    X->firstNode = v;
  else
    X->secondNode = addToSet(v,X->secondNode);
  return(X);
}

//node *makeNewNode(char *label, int i)
node *makeNewNode(int label, int i)
{
  return(makeNode(label,NULL,i));
}

//node *makeNode(char *label, edge *parentEdge, int index)
node *makeNode(int label, edge *parentEdge, int index)
{
  node *newNode;  /*points to new node added to the graph*/
  newNode = (node *) malloc(sizeof(node));
//  strncpy(newNode->label,label,NODE_LABEL_LENGTH);
  newNode->label = label;
  newNode->index = index;
  newNode->index2 = -1;
  newNode->parentEdge = parentEdge;
  newNode->leftEdge = NULL;
  newNode->middleEdge = NULL;
  newNode->rightEdge = NULL;
  /*all fields have been initialized*/
  return(newNode);
}

/*copyNode returns a copy of v which has all of the fields identical to those
of v, except the node pointer fields*/
node *copyNode(node *v)
{
  node *w;
  w = makeNode(v->label,NULL,v->index);
  w->index2 = v->index2;
  return(w);
}

edge *siblingEdge(edge *e)
{
  if(e == e->tail->leftEdge)
    return(e->tail->rightEdge);
  else
    return(e->tail->leftEdge);
}

edge *makeEdge(char *label, node *tail, node *head, double weight)
{
  edge *newEdge;
  newEdge = (edge *) malloc(sizeof(edge));
  strncpy(newEdge->label,label,EDGE_LABEL_LENGTH-1);
  newEdge->tail = tail;
  newEdge->head = head;
  newEdge->distance = weight;
  newEdge->totalweight = 0.0;
  return(newEdge);
}

tree *newTree()
{
  tree *T;
  T = (tree *) malloc(sizeof(tree));
  T->root = NULL;
  T->size = 0;
  T->weight = -1;
  return(T);
}

void updateSizes(edge *e, int direction)
{
  edge *f;
  switch(direction)
    {
    case UP:
      f = e->head->leftEdge;
      if (NULL != f)
	updateSizes(f,UP);
      f = e->head->rightEdge;
      if (NULL != f)
	updateSizes(f,UP);
      e->topsize++;
      break;
    case DOWN:
      f = siblingEdge(e);
      if (NULL != f)
	updateSizes(f,UP);
      f = e->tail->parentEdge;
      if (NULL != f)
	updateSizes(f,DOWN);
      e->bottomsize++;
      break;
    }
}

/*detrifurcate takes the (possibly trifurcated) input tree
  and reroots the tree to a leaf*/
/*assumes tree is only trifurcated at root*/
tree *detrifurcate(tree *T)
{
  node *v, *w;
  edge *e, *f;
  v = T->root;
  if(leaf(v))
    return(T);
  if (NULL != v->parentEdge)
    {
      error("root %d is poorly rooted.", v->label);
    }
  for(e = v->middleEdge, v->middleEdge = NULL; NULL != e; e = f )
    {
      w = e->head;
      v = e->tail;
      e->tail = w;
      e->head = v;
      f = w->leftEdge;
      v->parentEdge = e;
      w->leftEdge = e;
      w->parentEdge = NULL;
    }
  T->root = w;
  return(T);
}

void compareSets(tree *T, set *S)
{
  edge *e;
  node *v,*w;
  set *X;
  e = depthFirstTraverse(T,NULL);
  while (NULL != e)
    {
      v = e->head;
      for(X = S; NULL != X; X = X->secondNode)
	{
	  w = X->firstNode;
//	  if (0 == strcmp(v->label,w->label))
	  if (v->label == w->label)
	    {
	      v->index2 = w->index2;
	    w->index2 = -1;
	    break;
	    }
	}
      e = depthFirstTraverse(T,e);
    }
  v = T->root;
  for(X = S; NULL != X; X = X->secondNode)
    {
      w = X->firstNode;
//      if (0 == strcmp(v->label,w->label))
      if (v->label == w->label)
	{
	  v->index2 = w->index2;
	  w->index2 = -1;
	  break;
	}
    }
  if (-1 == v->index2)
    {
      error("leaf %d in tree not in distance matrix.", v->label);
    }
  e = depthFirstTraverse(T,NULL);
  while (NULL != e)
    {
      v = e->head;
      if ((leaf(v)) && (-1 == v->index2))
	{
	  error("leaf %d in tree not in distance matrix.", v->label);
	}
      e = depthFirstTraverse(T,e);
      }
  for(X = S; NULL != X; X = X->secondNode)
    if (X->firstNode->index2 > -1)
      {
	error("node %d in matrix but not a leaf in tree.", X->firstNode->label);
      }
  return;
}

void partitionSizes(tree *T)
{
  edge *e;
  e = depthFirstTraverse(T,NULL);
  while (NULL != e)
    {
      if (leaf(e->head))
	e->bottomsize = 1;
      else
	e->bottomsize = e->head->leftEdge->bottomsize
	  + e->head->rightEdge->bottomsize;
      e->topsize = (T->size + 2)/2 - e->bottomsize;
      e = depthFirstTraverse(T,e);
    }
}

/*************************************************************************

                           TRAVERSE FUNCTIONS

*************************************************************************/

edge *depthFirstTraverse(tree *T, edge *e)
     /*depthFirstTraverse returns the edge f which is least in T according
       to the depth-first order, but which is later than e in the search
       pattern.  If e is null, f is the least edge of T*/
{
  edge *f;
  if (NULL == e)
    {
      f = T->root->leftEdge;
      if (NULL != f)
	f = findBottomLeft(f);
      return(f);  /*this is the first edge of this search pattern*/
    }
  else /*e is non-null*/
    {
      if (e->tail->leftEdge == e)
	/*if e is a left-oriented edge, we skip the entire
	  tree cut below e, and find least edge*/
	f = moveRight(e);
      else  /*if e is a right-oriented edge, we have already looked at its
	      sibling and everything below e, so we move up*/
	f = e->tail->parentEdge;
    }
  return(f);
}

edge *findBottomLeft(edge *e)
     /*findBottomLeft searches by gottom down in the tree and to the left.*/
{
  edge *f;
  f = e;
  while (NULL != f->head->leftEdge)
    f = f->head->leftEdge;
  return(f);
}

edge *moveRight(edge *e)
{
  edge *f;
  f = e->tail->rightEdge; /*this step moves from a left-oriented edge
			    to a right-oriented edge*/
  if (NULL != f)
    f = findBottomLeft(f);
  return(f);
}

edge *topFirstTraverse(tree *T, edge *e)
     /*topFirstTraverse starts from the top of T, and from there moves stepwise
       down, left before right*/
     /*assumes tree has been detrifurcated*/
{
  edge *f;
  if (NULL == e)
    return(T->root->leftEdge); /*first Edge searched*/
  else if (!(leaf(e->head)))
    return(e->head->leftEdge); /*down and to the left is preferred*/
  else /*e->head is a leaf*/
    {
      f = moveUpRight(e);
      return(f);
    }
}

edge *moveUpRight(edge *e)
{
  edge *f;
  f = e;
  while ((NULL != f) && ( f->tail->leftEdge != f))
    f = f->tail->parentEdge;
  /*go up the tree until f is a leftEdge*/
  if (NULL == f)
    return(f); /*triggered at end of search*/
  else
    return(f->tail->rightEdge);
  /*and then go right*/
}

/*************************************************************************

                           FREE FUNCTIONS

*************************************************************************/

void freeMatrix(double **D, int size)
{
  int i;
  for(i=0;i<size;i++)
    free(D[i]);
  free(D);
}

void freeSet(set *S)
{
    if (NULL != S) {
	free(S->firstNode); /* added by EP 2014-03-04 */
	freeSet(S->secondNode);
    }
    free(S);
}

void freeTree(tree *T)
{
  node *v;
  v = T->root;
  if (NULL != v->leftEdge)
    freeSubTree(v->leftEdge);
  free(T->root);
  free(T);
}

void freeSubTree(edge *e)
{
  node *v;
  edge *e1, *e2;
  v = e->head;
  e1 = v->leftEdge;
  if (NULL != e1)
    freeSubTree(e1);
  e2 = v->rightEdge;
  if (NULL != e2)
    freeSubTree(e2);
  free(v);
  e->tail = NULL;
  e->head = NULL;
  free(e);
}
