/* me.h    2012-04-30 */

/* Copyright 2007-2008 Vincent Lefort, modified by Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>

#ifndef NONE
#define NONE 0
#endif
#ifndef UP
#define UP 1
#endif
#ifndef DOWN
#define DOWN 2
#endif
#ifndef LEFT
#define LEFT 3
#endif
#ifndef RIGHT
#define RIGHT 4
#endif
#ifndef SKEW
#define SKEW 5
#endif
#ifndef MAX_LABEL_LENGTH
#define MAX_LABEL_LENGTH 30
#endif
//#ifndef NODE_LABEL_LENGTH
//#define NODE_LABEL_LENGTH 30
//#endif
#ifndef EDGE_LABEL_LENGTH
#define EDGE_LABEL_LENGTH 30
#endif
#ifndef MAX_DIGITS
#define MAX_DIGITS 20
#endif
/* #ifndef INPUT_SIZE */
/* #define INPUT_SIZE 100 */
/* #endif */
#ifndef MAX_INPUT_SIZE
#define MAX_INPUT_SIZE 100000
#endif
#ifndef EPSILON
#define EPSILON 1.E-06
#endif
#ifndef ReadOpenParenthesis
#define ReadOpenParenthesis 0
#endif
#ifndef ReadSubTree
#define ReadSubTree 1
#endif
#ifndef ReadLabel
#define ReadLabel 2
#endif
#ifndef ReadWeight
#define ReadWeight 3
#endif
#ifndef AddEdge
#define AddEdge 4
#endif

#define XINDEX(i, j) n*(i - 1) - i*(i - 1)/2 + j - i - 1

typedef struct word
{
  char name[MAX_LABEL_LENGTH];
  struct word *suiv;
} WORD;

typedef struct pointers
{
  WORD *head;
  WORD *tail;
} POINTERS;

typedef struct node {
  int label; /* char label[NODE_LABEL_LENGTH]; */
  struct edge *parentEdge;
  struct edge *leftEdge;
  struct edge *middleEdge;
  struct edge *rightEdge;
  int index;
  int index2;
} node;

typedef struct edge {
  char label[EDGE_LABEL_LENGTH];
  struct node *tail; /*for edge (u,v), u is the tail, v is the head*/
  struct node *head;
  int bottomsize; /*number of nodes below edge */
  int topsize;    /*number of nodes above edge */
  double distance;
  double totalweight;
} edge;

typedef struct tree {
  char name[MAX_LABEL_LENGTH];
  struct node *root;
  int size;
  double weight;
} tree;

typedef struct set
{
  struct node *firstNode;
  struct set *secondNode;
} set;

void me_b(double *X, int *N, int *labels, int *nni, int *spr, int *tbr, int *edge1, int *edge2, double *el);
void me_o(double *X, int *N, int *labels, int *nni, int *edge1, int *edge2, double *el);
double **initDoubleMatrix(int d);
double **loadMatrix (double *X, int *labels, int n, set *S);
set *addToSet(node *v, set *X);
node *makeNewNode(int label, int i);
node *makeNode(int label, edge *parentEdge, int index);
node *copyNode(node *v);
edge *siblingEdge(edge *e);
edge *makeEdge(char *label, node *tail, node *head, double weight);
tree *newTree();
void updateSizes(edge *e, int direction);
tree *detrifurcate(tree *T);
void compareSets(tree *T, set *S);
void partitionSizes(tree *T);
edge *depthFirstTraverse(tree *T, edge *e);
edge *findBottomLeft(edge *e);
edge *moveRight(edge *e);
edge *topFirstTraverse(tree *T, edge *e);
edge *moveUpRight(edge *e);
void freeMatrix(double **D, int size);
void freeSet(set *S);
void freeTree(tree *T);
void freeSubTree(edge *e);
int leaf(node *v);
/* int whiteSpace(char c); */
/* node *decodeNewickSubtree(char *treeString, tree *T, int *uCount); */
/* tree *readNewickString (char *str, int numLeaves); */
void subtree2phylo(node *parent, int *edge1, int *edge2, double *el, int *ilab);
void tree2phylo(tree *T, int *edge1, int *edge2, double *el, int *ilab, int n);
