/* ape.c    2023-06-07 */

/* Copyright 2011-2023 Emmanuel Paradis, and 2007 R Development Core Team */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R_ext/Rdynload.h>
#include "ape.h"

int give_index(int i, int j, int n)
{
	if (i > j) return(DINDEX(j, i));
	else return(DINDEX(i, j));
}

/* From R-ext manual
   (not the same than in library/stats/src/nls.c) */
SEXP getListElement(SEXP list, char *str)
{
    SEXP res = R_NilValue, names = getAttrib(list, R_NamesSymbol);

    for (int i = 0; i < length(list); i++) {
	if (! strcmp(CHAR(STRING_ELT(names, i)), str)) {
	    res = VECTOR_ELT(list, i);
	    break;
	}
    }

    return res;
}

/* declare functions here to register them below */

void C_additive(double *dd, int* np, int* mp, double *ret);
void C_bionj(double *X, int *N, int *edge1, int *edge2, double *el);
void C_bionjs(double *D, int *N, int *edge1, int *edge2, double *edge_length, int* fsS);
void delta_plot(double *D, int *size, int *nbins, int *counts, double *deltabar);
SEXP dist_nodes(SEXP Ntips, SEXP Nnodes, SEXP edge, SEXP edge_length);
void C_ewLasso(double *D, int *N, int *e1, int *e2);
void mat_expo(double *P, int *nr);
void me_b(double *X, int *N, int *labels,
	  int *nni, int *spr, int *tbr,
	  int *edge1, int *edge2, double *el);
void me_o(double *X, int *N, int *labels, int *nni,
	  int *edge1, int *edge2, double *el);
void C_mvr(double *D, double* v,int *N, int *edge1, int *edge2, double *edge_length);
void C_mvrs(double *D, double* v, int *N, int *edge1, int *edge2, double *edge_length, int* fsS);
void neworder_pruningwise(int *ntip, int *nnode, int *edge1,
			  int *edge2, int *nedge, int *neworder);
SEXP C_nj(SEXP DIST, SEXP N);
void C_njs(double *D, int *N, int *edge1, int *edge2, double *edge_length, int *fsS);
void node_depth(int *ntip, int *e1, int *e2,
		int *nedge, double *xx, int *method);
void node_depth_edgelength(int *edge1, int *edge2, int *nedge,
			   double *edge_length, double *xx);
void node_height(int *edge1, int *edge2, int *nedge, double *yy);
void node_height_clado(int *ntip, int *edge1, int *edge2,
		       int *nedge, double *xx, double *yy);
void C_pic(int *ntip, int *edge1, int *edge2,
	   double *edge_len, double *phe, double *contr,
	   double *var_contr, int *var, int *scaled);
void C_rTraitCont(int *model, int *Nedge, int *edge1, int *edge2, double *el,
		  double *sigma, double *alpha, double *theta, double *x);
void C_treePop(int* splits, double* w,int* ncolp,int* np, int* ed1, int* ed2, double* edLen);
void C_triangMtd(double* d, int* np, int* ed1,int* ed2, double* edLen);
void C_triangMtds(double* d, int* np, int* ed1,int* ed2, double* edLen);
void C_ultrametric(double *dd, int* np, int* mp, double *ret);
void DNAbin2indelblock(unsigned char *x, int *n, int *s, int *y);
void trans_DNA2AA(unsigned char *x, int *s, unsigned char *res, int *code);

SEXP dist_dna(SEXP DNASEQ, SEXP MODEL, SEXP BASEFREQ, SEXP PAIRDEL, SEXP VARIANCE,
	      SEXP GAMMA, SEXP ALPHA);
SEXP GlobalDeletionDNA(SEXP DNASEQ);
SEXP C_where(SEXP DNASEQ, SEXP PAT);
SEXP rawStreamToDNAorAAbin(SEXP x, SEXP DNA);
SEXP seq_root2tip(SEXP edge, SEXP nbtip, SEXP nbnode);
/* SEXP treeBuildWithTokens(SEXP nwk); */
SEXP treeBuild(SEXP nwk);
/* SEXP cladoBuildWithTokens(SEXP nwk); */
SEXP cladoBuild(SEXP nwk);
SEXP bitsplits_multiPhylo(SEXP x, SEXP n, SEXP nr);
SEXP CountBipartitionsFromSplits(SEXP split, SEXP SPLIT);
SEXP _ape_prop_part2(SEXP trees, SEXP nTips);
SEXP _ape_bipartition2(SEXP orig, SEXP nTips);
SEXP _ape_reorderRcpp(SEXP orig, SEXP nTips, SEXP root, SEXP order);
SEXP writeDNAbinToFASTA(SEXP x, SEXP FILENAME, SEXP n, SEXP s, SEXP labels);
SEXP writeAAbinToFASTA(SEXP x, SEXP FILENAME, SEXP n, SEXP s, SEXP labels);
SEXP charVectorToDNAbinVector(SEXP x);
SEXP leading_trailing_gaps_to_N(SEXP DNASEQ);
SEXP SegSites(SEXP DNASEQ, SEXP STRICT);
SEXP BaseProportion(SEXP x);

static R_CMethodDef C_entries[] = {
    {"C_additive", (DL_FUNC) &C_additive, 4},
    {"C_bionj", (DL_FUNC) &C_bionj, 5},
    {"C_bionjs", (DL_FUNC) &C_bionjs, 6},
    {"delta_plot", (DL_FUNC) &delta_plot, 5},
    {"C_ewLasso", (DL_FUNC) &C_ewLasso, 4},
    {"mat_expo", (DL_FUNC) &mat_expo, 2},
    {"me_b", (DL_FUNC) &me_b, 9},
    {"me_o", (DL_FUNC) &me_o, 7},
    {"C_mvr", (DL_FUNC) &C_mvr, 6},
    {"C_mvrs", (DL_FUNC) &C_mvrs, 7},
    {"neworder_pruningwise", (DL_FUNC) &neworder_pruningwise, 6},
    {"C_njs", (DL_FUNC) &C_njs, 6},
    {"node_depth", (DL_FUNC) &node_depth, 6},
    {"node_depth_edgelength", (DL_FUNC) &node_depth_edgelength, 5},
    {"node_height", (DL_FUNC) &node_height, 4},
    {"node_height_clado", (DL_FUNC) &node_height_clado, 6},
    {"C_pic", (DL_FUNC) &C_pic, 9},
    {"C_rTraitCont", (DL_FUNC) &C_rTraitCont, 9},
    {"C_treePop", (DL_FUNC) &C_treePop, 7},
    {"C_triangMtd", (DL_FUNC) &C_triangMtd, 5},
    {"C_triangMtds", (DL_FUNC) &C_triangMtds, 5},
    {"C_ultrametric", (DL_FUNC) &C_ultrametric, 4},
    {"DNAbin2indelblock", (DL_FUNC) &DNAbin2indelblock, 4},
    {"trans_DNA2AA", (DL_FUNC) &trans_DNA2AA, 4},
    {NULL, NULL, 0}
};

static R_CallMethodDef Call_entries[] = {
    {"dist_dna", (DL_FUNC) &dist_dna, 7},
    {"dist_nodes", (DL_FUNC) &dist_nodes, 4},
    {"GlobalDeletionDNA", (DL_FUNC) &GlobalDeletionDNA, 1},
    {"rawStreamToDNAorAAbin", (DL_FUNC) &rawStreamToDNAorAAbin, 2},
    {"seq_root2tip", (DL_FUNC) &seq_root2tip, 3},
    /* {"treeBuildWithTokens", (DL_FUNC) &treeBuildWithTokens, 1}, */
    {"treeBuild", (DL_FUNC) &treeBuild, 1},
    /* {"cladoBuildWithTokens", (DL_FUNC) &cladoBuildWithTokens, 1}, */
    {"cladoBuild", (DL_FUNC) &cladoBuild, 1},
    {"bitsplits_multiPhylo", (DL_FUNC) &bitsplits_multiPhylo, 3},
    {"CountBipartitionsFromSplits", (DL_FUNC) &CountBipartitionsFromSplits, 2},
    {"BaseProportion", (DL_FUNC) &BaseProportion, 1},
    {"SegSites", (DL_FUNC) &SegSites, 2},
    {"C_where", (DL_FUNC) &C_where, 2},
    {"_ape_bipartition2", (DL_FUNC) &_ape_bipartition2, 2},
    {"_ape_prop_part2", (DL_FUNC) &_ape_prop_part2, 2},
    {"_ape_reorderRcpp", (DL_FUNC) &_ape_reorderRcpp, 4},
    {"writeDNAbinToFASTA", (DL_FUNC) &writeDNAbinToFASTA, 5},
    {"writeAAbinToFASTA", (DL_FUNC) &writeAAbinToFASTA, 5},
    {"charVectorToDNAbinVector", (DL_FUNC) &charVectorToDNAbinVector, 1},
    {"leading_trailing_gaps_to_N", (DL_FUNC) &leading_trailing_gaps_to_N, 1},
    {"C_nj", (DL_FUNC) &C_nj, 2},
    {NULL, NULL, 0}
};

void R_init_ape(DllInfo *info)
{
    R_registerRoutines(info, C_entries, Call_entries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
}
