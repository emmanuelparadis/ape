/* SPR.c    2013-09-26 */

/* Copyright 2009 Richard Desper */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "me.h"

/*functions from bNNI.c*/
void makeBMEAveragesTable(tree *T, double **D, double **A);
void assignBMEWeights(tree *T, double **A);

/*from me.c*/
edge *siblingEdge(edge *e);
void weighTree(tree *T);
void freeMatrix(double **D, int size);
edge *depthFirstTraverse(tree *T, edge *e);
double **initDoubleMatrix(int d);

/*from below*/
node *indexedNode(tree *T, int i);
edge *indexedEdge(tree *T, int i);
void assignSPRWeights(node *v, double **A, double ***swapWeights);
void SPRTopShift(tree *T, node *vmove, edge *esplit, int UpOrDown);
void assignDownWeightsUp(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights);
void assignDownWeightsSkew(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights);
void assignDownWeightsDown(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights);
void assignUpWeights(edge *etest, node *vtest, node *va, edge *back, node *cprve, double oldD_AB, double coeff, double **A, double ***swapWeights);

void zero3DMatrix(double ***X, int h, int l, int w)
{
	int i,j,k;
	for(i=0;i<h;i++)
		for(j=0;j<l;j++)
			for(k=0;k<w;k++)
				X[i][j][k] = 0.0;
}


void findTableMin(int *imin, int *jmin, int *kmin, int n, double ***X, double *min)
{
  int i,j,k;
  for(i=0;i<2;i++)
    for(j=0;j<n;j++)
      for(k=0;k<n;k++)
	{
	  if (X[i][j][k] < *min)
	    {
	      *min = X[i][j][k];
	      *imin = i;
	      *jmin = j;
	      *kmin = k;
	    }
	}
}


void SPR(tree *T, double **D, double **A, int *count)
{
  int i,j,k;
  node *v;
  /*FILE *treefile;*/
  edge *e,*f;
  /* char filename[MAX_LABEL_LENGTH];*/
  double ***swapWeights;
  double swapValue = 0.0;
  swapWeights = (double ***)malloc(2*sizeof(double **));
  makeBMEAveragesTable(T,D,A);
  assignBMEWeights(T,A);
  weighTree(T);
  /*if (verbose)
    printf("Before SPRs: tree length is %lf.\n",T->weight);*/
  for(i=0;i<2;i++)
    swapWeights[i] = initDoubleMatrix(T->size);
  do
    {
      swapValue=0.0;
      zero3DMatrix(swapWeights,2,T->size,T->size);
      i = j = k = 0;
      for(e=depthFirstTraverse(T,NULL);NULL!=e;e=depthFirstTraverse(T,e))
	assignSPRWeights(e->head,A,swapWeights);
      findTableMin(&i,&j,&k,T->size,swapWeights,&swapValue);
      swapValue = swapWeights[i][j][k];
      if (swapValue < -EPSILON)
	{
//	  if (verbose)
//	    printf("New tree weight should be %lf.\n",T->weight + 0.25*swapValue);
	  v = indexedNode(T,j);
	  f = indexedEdge(T,k);
//	  if (verbose)
//	    printf("Swapping tree below %s to split edge %s with head %s and tail %s\n",
//			   v->parentEdge->label,f->label,f->head->label,f->tail->label);
	  SPRTopShift(T,v,f,2-i);
	  makeBMEAveragesTable(T,D,A);
	  assignBMEWeights(T,A);
	  weighTree(T);
	  (*count)++;
	  /*sprintf(filename,"tree%d.new",*count);*/
//	  if (verbose)
//	    printf("After %d SPRs, tree weight is %lf.\n\n",*count,T->weight);
	  /*treefile = fopen(filename,"w");
	  NewickPrintTree(T,treefile);
	  fclose(treefile);*/
	  }
    } while (swapValue < -EPSILON);
  for(i=0;i<2;i++)
    freeMatrix(swapWeights[i],T->size);
  free(swapWeights);
  /*if (verbose)
    readOffTree(T);*/
}

/*assigns values to array swapWeights*/
/*swapWeights[0][j][k] will be the value of removing the tree below the edge whose head node has index j
and reattaching it to split the edge whose head has the index k*/
/*swapWeights[1][j][k] will be the value of removing the tree above the edge whose head node has index j
and reattaching it to split the edge whose head has the index k*/
void assignSPRWeights(node *vtest, double **A, double ***swapWeights)
{
  edge *etest, *left, *right, *sib, *par;
  etest = vtest->parentEdge;
  left = vtest->leftEdge;
  right = vtest->rightEdge;
  par = etest->tail->parentEdge;
  sib = siblingEdge(etest);
  if (NULL != par)
    assignDownWeightsUp(par,vtest,sib->head,NULL,NULL,0.0,1.0,A,swapWeights);
  if (NULL != sib)
    assignDownWeightsSkew(sib,vtest,sib->tail,NULL,NULL,0.0,1.0,A,swapWeights);
  /*assigns values for moving subtree rooted at vtest, starting with edge
    parental to tail of edge parental to vtest*/
  if (NULL != left)
    {
      assignUpWeights(left,vtest,right->head,NULL,NULL,0.0,1.0,A,swapWeights);
      assignUpWeights(right,vtest,left->head,NULL,NULL,0.0,1.0,A,swapWeights);
    }
}


/*recall NNI formula: change in tree length from AB|CD split to AC|BD split is
proportional to D_AC + D_BD - D_AB - D_CD*/
/*in our case B is the tree being moved (below vtest), A is the tree backwards below back, but
  with the vtest subtree removed, C is the sibling tree of back and D is the tree above etest*/
/*use va to denote the root of the sibling tree to B in the original tree*/
/*please excuse the multiple uses of the same letters: A,D, etc.*/
void assignDownWeightsUp(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights)
{
  edge *par, *sib, *skew;
  double D_AC, D_BD, D_AB, D_CD;
  par = etest->tail->parentEdge;
  skew = siblingEdge(etest);
  if (NULL == back) /*first recursive call*/
    {
      if (NULL == par)
	return;
      else /*start the process of assigning weights recursively*/
	{
	  assignDownWeightsUp(par,vtest,va,etest,va,A[va->index][vtest->index],0.5,A,swapWeights);
	  assignDownWeightsSkew(skew,vtest,va,etest,va,A[va->index][vtest->index],0.5,A,swapWeights);
	}
    }
  else /*second or later recursive call*/
    {
      sib = siblingEdge(back);
      D_BD = A[vtest->index][etest->head->index]; /*straightforward*/
      D_CD = A[sib->head->index][etest->head->index]; /*this one too*/
      D_AC = A[sib->head->index][back->head->index] + coeff*(A[sib->head->index][va->index]
							     - A[sib->head->index][vtest->index]);
      D_AB = 0.5*(oldD_AB + A[vtest->index][cprev->index]);
      swapWeights[0][vtest->index][etest->head->index] = swapWeights[0][vtest->index][back->head->index] + (D_AC + D_BD - D_AB - D_CD);
      if (NULL != par)
	{
	  assignDownWeightsUp(par,vtest,va,etest,sib->head,D_AB,0.5*coeff,A,swapWeights);
	  assignDownWeightsSkew(skew,vtest,va,etest,sib->head,D_AB,0.5*coeff,A,swapWeights);
	}
    }
}

void assignDownWeightsSkew(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights)
{
  /*same general idea as assignDownWeights, except needing to keep track of things a bit differently*/
  edge *par, *left, *right;
  /*par here = sib before
    left, right here = par, skew before*/
  double D_AB, D_CD, D_AC, D_BD;
  /*B is subtree being moved - below vtest
    A is subtree remaining fixed - below va, unioned with all trees already passed by B*/
  /*C is subtree being passed by B, in this case above par
    D is subtree below etest, fixed on other side*/
  par = etest->tail->parentEdge;
  left = etest->head->leftEdge;
  right = etest->head->rightEdge;
  if (NULL == back)
    {
      if (NULL == left)
	return;
      else
	{
	  assignDownWeightsDown(left,vtest,va,etest,etest->tail,A[vtest->index][etest->tail->index],0.5,A,swapWeights);
	  assignDownWeightsDown(right,vtest,va,etest,etest->tail,A[vtest->index][etest->tail->index],0.5,A,swapWeights);
	}
    }
  else
    {
      D_BD = A[vtest->index][etest->head->index];
      D_CD = A[par->head->index][etest->head->index];
      D_AC = A[back->head->index][par->head->index] + coeff*(A[va->index][par->head->index] - A[vtest->index][par->head->index]);
      D_AB = 0.5*(oldD_AB + A[vtest->index][cprev->index]);
      swapWeights[0][vtest->index][etest->head->index] = swapWeights[0][vtest->index][back->head->index] + (D_AC + D_BD - D_AB - D_CD);
      if (NULL != left)
	{
	  assignDownWeightsDown(left,vtest, va, etest, etest->tail, D_AB, 0.5*coeff, A, swapWeights);
	  assignDownWeightsDown(right,vtest, va, etest, etest->tail, D_AB, 0.5*coeff, A, swapWeights);
	}
    }
}

void assignDownWeightsDown(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A, double ***swapWeights)
{
  /*again the same general idea*/
  edge *sib, *left, *right;
  /*sib here = par in assignDownWeightsSkew
    rest is the same as assignDownWeightsSkew*/
  double D_AB, D_CD, D_AC, D_BD;
  /*B is below vtest, A is below va unioned with all trees already passed by B*/
  /*C is subtree being passed - below sib*/
  /*D is tree below etest*/
  sib = siblingEdge(etest);
  left = etest->head->leftEdge;
  right = etest->head->rightEdge;
  D_BD = A[vtest->index][etest->head->index];
  D_CD = A[sib->head->index][etest->head->index];
  D_AC = A[sib->head->index][back->head->index] + coeff*(A[sib->head->index][va->index] - A[sib->head->index][vtest->index]);
  D_AB = 0.5*(oldD_AB + A[vtest->index][cprev->index]);
  swapWeights[0][vtest->index][etest->head->index] = swapWeights[0][vtest->index][back->head->index] + ( D_AC + D_BD - D_AB - D_CD);
  if (NULL != left)
    {
      assignDownWeightsDown(left,vtest, va, etest, sib->head, D_AB, 0.5*coeff, A, swapWeights);
      assignDownWeightsDown(right,vtest, va, etest, sib->head, D_AB, 0.5*coeff, A, swapWeights);
    }
}

void assignUpWeights(edge *etest, node *vtest, node *va, edge *back, node *cprev, double oldD_AB, double coeff, double **A,
		     double ***swapWeights)
{
	/*SPR performed on tree above vtest...*/
	/*same idea as above, with appropriate selections of edges and nodes*/
  edge *sib, *left, *right;
  /*B is above vtest, A is other tree below vtest unioned with trees in path to vtest*/
  /*sib is tree C being passed by B*/
  /*D is tree below etest*/
  double D_AB, D_CD, D_AC, D_BD;
  // double thisWeight; deleted by EP, 2013-09-16, also below
  sib = siblingEdge(etest);
  left = etest->head->leftEdge;
  right = etest->head->rightEdge;
  if (NULL == back) /*first recursive call*/
    {
      if (NULL == left)
	return;
      else /*start the process of assigning weights recursively*/
	{
	  assignUpWeights(left,vtest,va,etest,va,A[va->index][vtest->index],0.5,A,swapWeights);
	  assignUpWeights(right,vtest,va,etest,va,A[va->index][vtest->index],0.5,A,swapWeights);
	}
    }
  else /*second or later recursive call*/
    {
      D_BD = A[vtest->index][etest->head->index];
      D_CD = A[sib->head->index][etest->head->index];
      D_AC = A[back->head->index][sib->head->index] + coeff*(A[va->index][sib->head->index] - A[vtest->index][sib->head->index]);
      D_AB = 0.5*(oldD_AB + A[vtest->index][cprev->index]);
      // thisWeight =  deleted by EP, 2013-09-16
      swapWeights[1][vtest->index][etest->head->index] = swapWeights[1][vtest->index][back->head->index] + (D_AC + D_BD - D_AB - D_CD);
      if (NULL != left)
	{
	  assignUpWeights(left,vtest, va, etest, sib->head, D_AB, 0.5*coeff, A, swapWeights);
	  assignUpWeights(right,vtest, va, etest, sib->head, D_AB, 0.5*coeff, A, swapWeights);
	}
    }
}

void pruneSubtree(edge *p, edge *u, edge *d)
/*starting with edge u above edges p, d*/
/*removes p, d from tree, u connects to d->head to compensate*/
{
  p->tail->parentEdge = NULL;/*remove p subtree*/
  u->head = d->head;
  d->head->parentEdge = u;	/*u connected to d->head*/
  d->head = NULL; /*d removed from tree*/
}

void SPRsplitEdge(edge *e, edge *p, edge *d)
/*splits edge e to make it parental to p,d.  d is parental to what
  previously was below e*/
{
  d->head = e->head;
  e->head = p->tail;
  p->tail->parentEdge = e;
  d->head->parentEdge = d;
}


/*topological shift function*/
/*removes subtree rooted at v and re-inserts to spilt e*/
void SPRDownShift(tree *T, node *v, edge *e)
{
  edge *vup, *vdown, *vpar;
  vpar = v->parentEdge;
  vdown = siblingEdge(vpar);
  vup = vpar->tail->parentEdge;
  /*topological shift*/
  pruneSubtree(vpar,vup,vdown);
  /*removes v subtree and vdown, extends vup*/
  SPRsplitEdge(e,vpar,vdown); /*splits e to make e sibling edge to vpar,
				both below vup*/
}

void SPRUpShift(tree *T, node *vmove, edge *esplit)
/*an inelegant iterative version*/
{
  edge *f;
  edge **EPath;
  edge **sib;
  node **v;
  int i;
  int pathLength;
  for(f=esplit->tail->parentEdge,pathLength=1;f->tail != vmove;f=f->tail->parentEdge)
    pathLength++;
  /*count number of edges to vmove*/
  /*note that pathLength > 0 will hold*/

  EPath = (edge **)malloc(pathLength*sizeof(edge *));
  v = (node **)malloc(pathLength*sizeof(edge *));
  sib = (edge **)malloc((pathLength+1)*sizeof(edge *));
  /*there are pathLength + 1 side trees, one at the head and tail of each edge in the path*/

  sib[pathLength] = siblingEdge(esplit);
  i = pathLength;
  f = esplit->tail->parentEdge;
  while (i > 0)
    {
      i--;
      EPath[i] = f;
      sib[i] = siblingEdge(f);
      v[i] = f->head;
      f = f->tail->parentEdge;
    }
  /*indexed so head of Epath is v[i], tail is v[i-1] and sibling edge is sib[i]*/
  /*need to assign head, tail of each edge in path
    as well as have proper values for the left and right fields*/

  if (esplit == esplit->tail->leftEdge)
    {
      vmove->leftEdge = esplit;
      vmove->rightEdge = EPath[pathLength-1];
    }
  else
    {
      vmove->rightEdge = esplit;
      vmove->leftEdge = EPath[pathLength-1];
    }
  esplit->tail = vmove;
  /*espilt->head remains unchanged*/
  /*vmove has the proper fields for left, right, and parentEdge*/

  for(i=0;i<(pathLength-1);i++)
    EPath[i]->tail = v[i+1];

  /*this bit flips the orientation along the path
    the tail of Epath[i] is now v[i+1] instead of v[i-1]*/

  EPath[pathLength-1]->tail = vmove;

  for(i=1;i<pathLength;i++)
    {
      if (sib[i+1] == v[i]->leftEdge)
	v[i]->rightEdge = EPath[i-1];
      else
	v[i]->leftEdge = EPath[i-1];
    }
  if (sib[1] == v[0]->leftEdge)
    v[0]->rightEdge = sib[0];
  else
    v[0]->leftEdge = sib[0];
  sib[0]->tail = v[0];
  free(EPath);
  free(v);
  free(sib);
}


void SPRTopShift(tree *T, node *vmove, edge *esplit, int UpOrDown)
{
  if (DOWN == UpOrDown)
    SPRDownShift(T,vmove,esplit);
  else
    SPRUpShift(T,vmove,esplit);
}

node *indexedNode(tree *T, int i)
{
  edge *e;
  for(e = depthFirstTraverse(T,NULL);NULL!=e;e=depthFirstTraverse(T,e))
    if (i == e->head->index)
      return(e->head);
  if (i == T->root->index)
    return(T->root);
  return(NULL);
}

edge *indexedEdge(tree *T, int i)
{
  edge *e;
  for(e = depthFirstTraverse(T,NULL);NULL!=e;e=depthFirstTraverse(T,e))
    if (i == e->head->index)
      return(e);
  return(NULL);
}
