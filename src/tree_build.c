/* tree_build.c    2020-02-12 */

/* Copyright 2008-2020 Emmanuel Paradis, 2017 Klaus Schliep */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R.h>
#include <Rinternals.h>

static int str2int(char *x, int n)
{
	int i, k = 1, ans = 0;

	for (i = n - 1; i >= 0; i--, k *= 10)
		ans += ((int)x[i] - 48) * k;

	return ans;
}

void extract_portion_Newick(const char *x, int a, int b, char *y)
{
	int i, j;

	for (i = a, j = 0; i <= b; i++, j++) y[j] = x[i];

	y[j] = '\0';
}

void decode_terminal_edge_token(const char *x, int a, int b, int *node, double *w)
{
	int co = a;
	char *endstr, str[100];

	while (x[co] != ':' && co <= b) co++;

	extract_portion_Newick(x, a, co - 1, str);
	*node = str2int(str, co - a);
        if (co < b) {
	    extract_portion_Newick(x, co + 1, b, str);
	    *w = R_strtod(str, &endstr);
        } else *w = NAN;
}

void decode_internal_edge(const char *x, int a, int b, char *lab, double *w)
{
	int co = a;
	char *endstr, str[100];

	while (x[co] != ':' && co <= b) co++;

	if (a == co) lab[0] = '\0'; /* if no node label */
	else extract_portion_Newick(x, a, co - 1, lab);
        if (co < b) {
	    extract_portion_Newick(x, co + 1, b, str);
	    *w = R_strtod(str, &endstr);
        } else *w = NAN;
}

void decode_terminal_edge_token_clado(const char *x, int a, int b, int *node)
{
	char str[100];  // *endstr,

	extract_portion_Newick(x, a, b, str);
	*node = str2int(str, b + 1 - a);
}

void decode_internal_edge_clado(const char *x, int a, int b, char *lab)
{
//	char *endstr, str[100];
	if (a > b) lab[0] = '\0'; /* if no node label */
	else extract_portion_Newick(x, a, b, lab);
}

void decode_terminal_edge(const char *x, int a, int b, char *tip, double *w)
{
	int co = a;
	char *endstr, str[100];

	while (x[co] != ':' && co <= b) co++;

	extract_portion_Newick(x, a, co - 1, tip);
	if (co < b) {
           extract_portion_Newick(x, co + 1, b, str);
           *w = R_strtod(str, &endstr);
        } else *w = NAN;
}

void decode_terminal_edge_clado(const char *x, int a, int b, char *tip)
{
	extract_portion_Newick(x, a, b, tip);
}

#define ADD_INTERNAL_EDGE            \
    e[j] = curnode;                  \
    e[j + nedge] = curnode = ++node; \
    stack_internal[k++] = j;         \
    j++

#define ADD_TERMINAL_EDGE                                        \
    e[j] = curnode;                                              \
    decode_terminal_edge_token(x, pr + 1, ps - 1, &tmpi, &tmpd); \
    e[j + nedge] = tmpi;                                         \
    el[j] = tmpd;                                                \
    j++

#define GO_DOWN                                                  \
    decode_internal_edge(x, ps + 1, pt - 1, lab, &tmpd);         \
    SET_STRING_ELT(node_label, curnode - 1 - ntip, mkChar(lab)); \
    l = stack_internal[--k];					 \
    el[l] = tmpd;                                                \
    curnode = e[l]

#define ADD_TERMINAL_EDGE_CLADO                                  \
    e[j] = curnode;                                              \
    decode_terminal_edge_token_clado(x, pr + 1, ps - 1, &tmpi);  \
    e[j + nedge] = tmpi;                                         \
    j++

#define GO_DOWN_CLADO                                            \
    decode_internal_edge_clado(x, ps + 1, pt - 1, lab);          \
    SET_STRING_ELT(node_label, curnode - 1 - ntip, mkChar(lab)); \
    l = stack_internal[--k];					 \
    curnode = e[l]


#define ADD_TERMINAL_EDGE_TIPLABEL                       \
    e[j] = curnode;                                      \
    decode_terminal_edge(x, pr + 1, ps - 1, tip, &tmpd); \
    SET_STRING_ELT(tip_label, curtip-1, mkChar(tip));    \
    e[j + nedge] = curtip;                               \
    el[j] = tmpd;                                        \
    curtip++;                                            \
    j++

#define ADD_TERMINAL_EDGE_TIPLABEL_CLADO                \
    e[j] = curnode;                                     \
    decode_terminal_edge_clado(x, pr + 1, ps - 1, tip); \
    SET_STRING_ELT(tip_label, curtip-1, mkChar(tip));   \
    e[j + nedge] = curtip;                              \
    curtip++;                                           \
    j++

#define INITIALIZE_SKELETON			 \
    PROTECT(nwk = coerceVector(nwk, STRSXP));	 \
    x = CHAR(STRING_ELT(nwk, 0));		 \
    n = strlen(x);				 \
    skeleton = (int *)R_alloc(n, sizeof(int *)); \
    for (i = 0; i < n; i++) {			 \
	if (x[i] == '(') {			 \
	    skeleton[nsk] = i;			 \
	    nsk++;				 \
	    nleft++;				 \
	    continue;				 \
	}					 \
	if (x[i] == ',') {			 \
	    skeleton[nsk] = i;			 \
	    nsk++;				 \
	    ntip++;				 \
	    continue;				 \
	}					 \
	if (x[i] == ')') {			 \
	    skeleton[nsk] = i;			 \
	    nsk++;				 \
	    nright++;				 \
	    nnode++;				 \
	}					 \
    }						 \
    if (nleft != nright) error("numbers of left and right parentheses in Newick string not equal\n"); \
    nedge = ntip + nnode - 1

/*
   NOTE: the four functions below use the same algorithm to build a
   "phylo" object from a Newick string (with/without edge lengths
   and/or with/without tokens). Only the first one is commented.

*/

SEXP treeBuildWithTokens(SEXP nwk)
{
	const char *x;
	int n, i, ntip = 1, nleft = 0, nright = 0, nnode = 0, nedge, *e, curnode, node, j, *skeleton, nsk = 0, ps, pr, pt, tmpi, l, k, stack_internal[10000];
	double *el, tmpd;
	char lab[512];
	SEXP edge, edge_length, Nnode, node_label, phy;

	/* first pass on the Newick string to localize parentheses and commas */
	INITIALIZE_SKELETON;

	PROTECT(Nnode = allocVector(INTSXP, 1));
	PROTECT(edge = allocVector(INTSXP, nedge*2));
	PROTECT(edge_length = allocVector(REALSXP, nedge));
	PROTECT(node_label = allocVector(STRSXP, nnode));
	INTEGER(Nnode)[0] = nnode;

	e = INTEGER(edge);
	el = REAL(edge_length);

	curnode = node = ntip + 1;
	k = j = 0;
/* j: index of the current position in the edge matrix */
/* k: index of the current position in stack_internal */
/* stack_internal is a simple array storing the indices of the
   successive internal edges from the root; it's a stack so it is
   incremented every time an internal edge is added, and decremented
   every GO_DOWN step. This makes easy to find the index of the
   subtending edge. */

/* second pass on the Newick string to build the "phylo" object elements */
	for (i = 1; i < nsk - 1; i++) {
		ps = skeleton[i];
		if (x[ps] == '(') {
			ADD_INTERNAL_EDGE;
			continue;
		}
		pr = skeleton[i - 1];
		if (x[ps] == ',') {
			if (x[pr] != ')') {
				/* !!! accolades indispensables !!! */
				ADD_TERMINAL_EDGE;
			}
			continue;
		}
		if (x[ps] == ')') {
			pt = skeleton[i + 1]; // <- utile ???
			if (x[pr] == ',') {
				ADD_TERMINAL_EDGE;
 				GO_DOWN;
				continue;
			}
			/* added by Klaus to allow singleton nodes (2017-05-28): */
			if (x[pr] == '(') {
			    ADD_TERMINAL_EDGE;
			    GO_DOWN;
			    continue;
			} /* end */
			if (x[pr] == ')') {
				GO_DOWN;
			}
		}
	}

	pr = skeleton[nsk - 2];
	ps = skeleton[nsk - 1];
/* is the last edge terminal? */
	if (x[pr] == ',' && x[ps] == ')') {
		ADD_TERMINAL_EDGE;
	}

/* is there a root edge and/or root label? */
	if (ps < n - 2) {
		i = ps + 1;
		while (i < n - 2 && x[i] != ':') i++;
		if (i < n - 2) {
			PROTECT(phy = allocVector(VECSXP, 5));
			SEXP root_edge;
			decode_internal_edge(x, ps + 1, n - 2, lab, &tmpd);
			PROTECT(root_edge = allocVector(REALSXP, 1));
			REAL(root_edge)[0] = tmpd;
			SET_VECTOR_ELT(phy, 4, root_edge);
			UNPROTECT(1);
			SET_STRING_ELT(node_label, 0, mkChar(lab));
		} else {
			extract_portion_Newick(x, ps + 1, n - 2, lab);
			SET_STRING_ELT(node_label, 0, mkChar(lab));
			PROTECT(phy = allocVector(VECSXP, 4));
		}
	} else PROTECT(phy = allocVector(VECSXP, 4));

	SET_VECTOR_ELT(phy, 0, edge);
	SET_VECTOR_ELT(phy, 1, edge_length);
	SET_VECTOR_ELT(phy, 2, Nnode);
	SET_VECTOR_ELT(phy, 3, node_label);

	UNPROTECT(6);
	return phy;
}

SEXP cladoBuildWithTokens(SEXP nwk)
{
	const char *x;
	int n, i, ntip = 1, nleft = 0, nright = 0, nnode = 0, nedge, *e, curnode, node, j, *skeleton, nsk = 0, ps, pr, pt, tmpi, l, k, stack_internal[10000];
	char lab[512];
	SEXP edge, Nnode, node_label, phy;

	INITIALIZE_SKELETON;

	PROTECT(Nnode = allocVector(INTSXP, 1));
	PROTECT(edge = allocVector(INTSXP, nedge*2));
	PROTECT(node_label = allocVector(STRSXP, nnode));
	INTEGER(Nnode)[0] = nnode;

	e = INTEGER(edge);

	curnode = node = ntip + 1;
	k = j = 0;

	for (i = 1; i < nsk - 1; i++) {
		ps = skeleton[i];
		if (x[ps] == '(') {
			ADD_INTERNAL_EDGE;
			continue;
		}
		pr = skeleton[i - 1];
		if (x[ps] == ',') {
			if (x[pr] != ')') {
				ADD_TERMINAL_EDGE_CLADO;
			}
			continue;
		}
		if (x[ps] == ')') {
			pt = skeleton[i + 1];
			if (x[pr] == ',') {
				ADD_TERMINAL_EDGE_CLADO;
 				GO_DOWN_CLADO;
				continue;
			}
			if (x[pr] == '(') {
			    ADD_TERMINAL_EDGE_CLADO;
			    GO_DOWN_CLADO;
			    continue;
			}
			if (x[pr] == ')') {
				GO_DOWN_CLADO;
			}
		}
	}

	pr = skeleton[nsk - 2];
	ps = skeleton[nsk - 1];
	if (x[pr] == ',' && x[ps] == ')') {
		ADD_TERMINAL_EDGE_CLADO;
	}
	if (ps < n - 2) {
	    extract_portion_Newick(x, ps + 1, n - 2, lab);
	    SET_STRING_ELT(node_label, 0, mkChar(lab));
	    PROTECT(phy = allocVector(VECSXP, 3));
	} else PROTECT(phy = allocVector(VECSXP, 3));

	SET_VECTOR_ELT(phy, 0, edge);
	SET_VECTOR_ELT(phy, 1, Nnode);
	SET_VECTOR_ELT(phy, 2, node_label);

	UNPROTECT(5);
	return phy;
}

SEXP treeBuild(SEXP nwk)
{
	const char *x;
	int n, i, ntip = 1, nleft = 0, nright = 0, nnode = 0, nedge, *e, curnode, node, j, *skeleton, nsk = 0, ps, pr, pt, l, k, stack_internal[10000], curtip = 1;
	double *el, tmpd;
	char lab[512], tip[512];
	SEXP edge, edge_length, Nnode, node_label, tip_label, phy;

	INITIALIZE_SKELETON;

	PROTECT(Nnode = allocVector(INTSXP, 1));
	PROTECT(edge = allocVector(INTSXP, nedge*2));
	PROTECT(edge_length = allocVector(REALSXP, nedge));
	PROTECT(node_label = allocVector(STRSXP, nnode));
        PROTECT(tip_label = allocVector(STRSXP, ntip));
	INTEGER(Nnode)[0] = nnode;

	e = INTEGER(edge);
	el = REAL(edge_length);

	curnode = node = ntip + 1;
	k = j = 0;

	for (i = 1; i < nsk - 1; i++) {
		ps = skeleton[i];
		if (x[ps] == '(') {
			ADD_INTERNAL_EDGE;
			continue;
		}
		pr = skeleton[i - 1];
		if (x[ps] == ',') {
			if (x[pr] != ')') {
				ADD_TERMINAL_EDGE_TIPLABEL;
			}
			continue;
		}
		if (x[ps] == ')') {
			pt = skeleton[i + 1];
			if (x[pr] == ',') {
				ADD_TERMINAL_EDGE_TIPLABEL;
 				GO_DOWN;
				continue;
			}
			if (x[pr] == '(') {
			    ADD_TERMINAL_EDGE_TIPLABEL;
			    GO_DOWN;
			    continue;
			}
			if (x[pr] == ')') {
				GO_DOWN;
			}
		}
	}

	pr = skeleton[nsk - 2];
	ps = skeleton[nsk - 1];
	if (x[pr] == ',' && x[ps] == ')') {
		ADD_TERMINAL_EDGE_TIPLABEL;
	}

	if (ps < n - 2) {
		i = ps + 1;
		while (i < n - 2 && x[i] != ':') i++;
		if (i < n - 2) {
			PROTECT(phy = allocVector(VECSXP, 6));
			SEXP root_edge;
			decode_internal_edge(x, ps + 1, n - 2, lab, &tmpd);
			PROTECT(root_edge = allocVector(REALSXP, 1));
			REAL(root_edge)[0] = tmpd;
			SET_VECTOR_ELT(phy, 5, root_edge);
			UNPROTECT(1);
			SET_STRING_ELT(node_label, 0, mkChar(lab));
		} else {
			extract_portion_Newick(x, ps + 1, n - 2, lab);
			SET_STRING_ELT(node_label, 0, mkChar(lab));
			PROTECT(phy = allocVector(VECSXP, 5));
		}
	} else PROTECT(phy = allocVector(VECSXP, 5));

	SET_VECTOR_ELT(phy, 0, edge);
	SET_VECTOR_ELT(phy, 1, edge_length);
	SET_VECTOR_ELT(phy, 2, Nnode);
	SET_VECTOR_ELT(phy, 3, node_label);
        SET_VECTOR_ELT(phy, 4, tip_label);
	UNPROTECT(7);
	return phy;
}

SEXP cladoBuild(SEXP nwk)
{
	const char *x;
	int n, i, ntip = 1, nleft = 0, nright = 0, nnode = 0, nedge, *e, curnode, node, j, *skeleton, nsk = 0, ps, pr, pt, l, k, stack_internal[10000], curtip = 1;
	char lab[512], tip[512];
	SEXP edge, Nnode, node_label, tip_label, phy;

	INITIALIZE_SKELETON;

	PROTECT(Nnode = allocVector(INTSXP, 1));
	PROTECT(edge = allocVector(INTSXP, nedge*2));
	PROTECT(node_label = allocVector(STRSXP, nnode));
        PROTECT(tip_label = allocVector(STRSXP, ntip));
	INTEGER(Nnode)[0] = nnode;

	e = INTEGER(edge);

	curnode = node = ntip + 1;
	k = j = 0;

	for (i = 1; i < nsk - 1; i++) {
		ps = skeleton[i];
		if (x[ps] == '(') {
			ADD_INTERNAL_EDGE;
			continue;
		}
		pr = skeleton[i - 1];
		if (x[ps] == ',') {
			if (x[pr] != ')') {
				ADD_TERMINAL_EDGE_TIPLABEL_CLADO;
			}
			continue;
		}
		if (x[ps] == ')') {
			pt = skeleton[i + 1];
			if (x[pr] == ',') {
				ADD_TERMINAL_EDGE_TIPLABEL_CLADO;
 				GO_DOWN_CLADO;
				continue;
			}
			if (x[pr] == '(') {
			    ADD_TERMINAL_EDGE_TIPLABEL_CLADO;
			    GO_DOWN_CLADO;
			    continue;
			}
			if (x[pr] == ')') {
			    GO_DOWN_CLADO;
			}
		}
	}

	pr = skeleton[nsk - 2];
	ps = skeleton[nsk - 1];
	if (x[pr] == ',' && x[ps] == ')') {
		ADD_TERMINAL_EDGE_TIPLABEL_CLADO;
	}

	if (ps < n - 2) {
			extract_portion_Newick(x, ps + 1, n - 2, lab);
			SET_STRING_ELT(node_label, 0, mkChar(lab));
			PROTECT(phy = allocVector(VECSXP, 4));
	} else PROTECT(phy = allocVector(VECSXP, 4));

	SET_VECTOR_ELT(phy, 0, edge);
	SET_VECTOR_ELT(phy, 1, Nnode);
	SET_VECTOR_ELT(phy, 2, node_label);
        SET_VECTOR_ELT(phy, 3, tip_label);
	UNPROTECT(6);
	return phy;
}

#undef ADD_INTERNAL_EDGE
#undef ADD_TERMINAL_EDGE
#undef ADD_TERMINAL_EDGE_CLADO
#undef ADD_TERMINAL_EDGE_TIPLABEL
#undef ADD_TERMINAL_EDGE_TIPLABEL_CLADO
#undef GO_DOWN
#undef GO_DOWN_CLADO
