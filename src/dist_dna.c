/* dist_dna.c       2020-08-18 */

/* Copyright 2005-2020 Emmanuel Paradis */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <R_ext/Lapack.h>
#include "ape.h"

/* from R: print(log(4), d = 22) */
#define LN4 1.386294361119890572454

/* returns 8 if the base is known surely, 0 otherwise */
#define KnownBase(a) (a & 8)

/* returns 1 if the base is adenine surely, 0 otherwise */
#define IsAdenine(a) (a == 136)

/* returns 1 if the base is guanine surely, 0 otherwise */
#define IsGuanine(a) (a == 72)

/* returns 1 if the base is cytosine surely, 0 otherwise */
#define IsCytosine(a) (a == 40)

/* returns 1 if the base is thymine surely, 0 otherwise */
#define IsThymine(a) (a == 24)

/* returns 1 if the base is a purine surely, 0 otherwise */
#define IsPurine(a) (a > 63)

/* returns 1 if the base is a pyrimidine surely, 0 otherwise */
#define IsPyrimidine(a) (a < 64)

/* returns 1 if both bases are different surely, 0 otherwise */
#define DifferentBase(a, b) ((a & b) < 16)

/* returns 1 if both bases are the same surely, 0 otherwise */
#define SameBase(a, b) (KnownBase(a) && a == b)

/* computes directly the determinant of a 4x4 matrix */
double detFourByFour(double *x)
{
    double det, a33a44, a34a43, a34a42, a32a44, a32a43, a33a42, a34a41, a31a44, a31a43, a33a41, a31a42, a32a41;

    a33a44 = x[10]*x[15]; a34a43 = x[14]*x[11];
    a34a42 = x[14]*x[7];  a32a44 = x[6]*x[15];
    a32a43 = x[6]*x[11];  a33a42 = x[10]*x[7];
    a34a41 = x[14]*x[3];  a31a44 = x[2]*x[15];
    a31a43 = x[2]*x[11];  a33a41 = x[10]*x[3];
    a31a42 = x[2]*x[7];   a32a41 = x[6]*x[3];

    det = x[0]*x[5]*(a33a44 - a34a43) + x[0]*x[9]*(a34a42 - a32a44) +
      x[0]*x[13]*(a32a43 - a33a42) + x[4]*x[9]*(a31a44 - a34a41) +
      x[4]*x[13]*(a33a41 - a31a43) + x[4]*x[1]*(a34a43 - a33a44) +
      x[8]*x[13]*(a31a42 - a32a41) + x[8]*x[1]*(a32a44 - a34a42) +
      x[8]*x[5]*(a34a41 - a31a44) + x[12]*x[1]*(a33a42 - a32a43) +
      x[12]*x[5]*(a31a43 - a33a41) + x[12]*x[9]*(a32a41 - a31a42);

    return det;
}

#define CHECK_PAIRWISE_DELETION\
    if (KnownBase(x[s1]) && KnownBase(x[s2])) L++;\
    else continue;

#define COUNT_TS_TV\
    if (SameBase(x[s1], x[s2])) continue;\
    Nd++;\
    if (IsPurine(x[s1]) && IsPurine(x[s2])) {\
        Ns++;\
        continue;\
    }\
    if (IsPyrimidine(x[s1]) && IsPyrimidine(x[s2])) Ns++;

void distDNA_indel(unsigned char *x, int n, int s, double *d)
{
	int i1, i2, s1, s2, target, N;

	target = 0;
	for (i1 = 1; i1 < n; i1++) {
		for (i2 = i1 + 1; i2 <= n; i2++) {
			N = 0;

			for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n * (s - 1); s1 += n, s2 += n)
				if ((x[s1] ^ x[s2]) & 4) N++;

			d[target] = ((double) N);
			target++;
		}
	}
}

void DNAbin2indelblock(unsigned char *x, int *n, int *s, int *y)
{
    int i, j, k, pos, ngap, indel = 0;

    for (i = 0; i < *n; i++) {
	j = i;
	k = 0;
	while (k < *s) {
	    if (x[j] == 4) {
		if (!indel) {
		    pos = j;
		    indel = 1;
		    ngap = 1;
		} else ngap++;
	    } else {
		if (indel) {
		    y[pos] = ngap;
		    indel = 0;
		}
	    }
	    j += *n;
	    k++;
	}
	if (indel) {
	    y[pos] = ngap;
	    indel = 0;
	}
    }
}

void distDNA_indelblock(unsigned char *x, int n, int s, double *d)
{
    int *y, i1, i2, s1, s2, target, Nd;

    y = (int*)R_alloc(n * s, sizeof(int));
    memset(y, 0, n * s * sizeof(int));
    DNAbin2indelblock(x, &n, &s, y);

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
	for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n * (s - 1); s1 += n, s2 += n)
		if (y[s1] != y[s2]) Nd++;
	    d[target] = ((double) Nd);
	    target++;
	}
    }
}

void distDNA_TsTv(unsigned char *x, int n, int s, double *d, int Ts, int pairdel)
{
	int i1, i2, s1, s2, target, Nd, Ns;

	target = 0;
	for (i1 = 1; i1 < n; i1++) {
		for (i2 = i1 + 1; i2 <= n; i2++) {
			Nd = Ns = 0;
			for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
				if (pairdel && !(KnownBase(x[s1]) && KnownBase(x[s2]))) continue;
				COUNT_TS_TV
			}
			if (Ts) d[target] = ((double) Ns); /* output number of transitions */
			else d[target] = ((double) Nd - Ns); /* output number of transversions */
			target++;
		}
	}
}

void distDNA_raw(unsigned char *x, int n, int s, double *d, int scaled)
{
	int i1, i2, s1, s2, target, Nd;

	target = 0;
	for (i1 = 1; i1 < n; i1++) {
		for (i2 = i1 + 1; i2 <= n; i2++) {
			Nd = 0;
			for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n * (s - 1); s1 += n, s2 += n)
				if (DifferentBase(x[s1], x[s2])) Nd++;
			if (scaled) d[target] = ((double) Nd / s);
			else d[target] = ((double) Nd);
			target++;
		}
	}
}

void distDNA_raw_pairdel(unsigned char *x, int n, int s, double *d, int scaled)
{
    int i1, i2, s1, s2, target, Nd, L;

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n * (s - 1); s1 += n, s2 += n) {
                CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    if (scaled) d[target] = ((double) Nd/L);
	    else d[target] = ((double) Nd);
	    target++;
	}
    }
}

#define COMPUTE_DIST_JC69\
    p = ((double) Nd/L);\
    if (gamma)\
      d[target] = 0.75 * alpha * (pow(1 - 4*p/3, -1/alpha) - 1);\
    else d[target] = -0.75 * log(1 - 4*p/3);\
    if (variance) {\
        if (gamma) var[target] = p*(1 - p)/(pow(1 - 4*p/3, -2/(alpha + 1)) * L);\
	else var[target] = p*(1 - p)/(pow(1 - 4*p/3, 2)*L);\
    }

void distDNA_JC69(unsigned char *x, int n, int s, double *d,
		  int variance, double *var, int gamma, double alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p;

    L = s;

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
  	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n * (s - 1); s1 += n, s2 += n)
	      if (DifferentBase(x[s1], x[s2])) Nd++;
	    COMPUTE_DIST_JC69
	    target++;
	}
    }
}

void distDNA_JC69_pairdel(unsigned char *x, int n, int s, double *d,
			  int variance, double *var, int gamma, double alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p;

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
  	    Nd = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1+= n, s2 += n) {
	        CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    COMPUTE_DIST_JC69
	    target++;
	}
    }
}

#define COMPUTE_DIST_K80\
    P = ((double) Ns/L);\
    Q = ((double) (Nd - Ns)/L);\
    a1 = 1 - 2*P - Q;\
    a2 = 1 - 2*Q;\
    if (gamma) {\
        b = -1 / alpha;\
    	d[target] = alpha * (pow(a1, b) + 0.5*pow(a2, b) - 1.5)/2;\
    }\
    else d[target] = -0.5 * log(a1 * sqrt(a2));\
    if (variance) {\
        if (gamma) {\
    	    b = -(1 / alpha + 1);\
    	    c1 = pow(a1, b);\
    	    c2 = pow(a2, b);\
    	    c3 = (c1 + c2)/2;\
    	} else {\
    	  c1 = 1/a1;\
    	  c2 = 1/a2;\
    	  c3 = (c1 + c2)/2;\
    	}\
    	var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;\
    }

void distDNA_K80(unsigned char *x, int n, int s, double *d,
		 int variance, double *var, int gamma, double alpha)
{
    int i1, i2, s1, s2, target, Nd, Ns, L;
    double P, Q, a1, a2, b, c1, c2, c3;

    L = s;

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1+= n, s2 += n) {
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_K80
	    target++;
	}
    }
}

void distDNA_K80_pairdel(unsigned char *x, int n, int s, double *d,
			 int variance, double *var, int gamma, double alpha)
{
    int i1, i2, s1, s2, target, Nd, Ns, L;
    double P, Q, a1, a2, b, c1, c2, c3;

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        CHECK_PAIRWISE_DELETION
		COUNT_TS_TV
	    }
	    COMPUTE_DIST_K80
	    target++;
	}
    }
}

#define COMPUTE_DIST_F81\
    p = ((double) Nd/L);\
    if (gamma) d[target] = E * alpha * (pow(1 - p/E, -1/ alpha) - 1);\
    else d[target] = -E*log(1 - p/E);\
    if (variance) {\
	if (gamma) var[target] = p*(1 - p)/(pow(1 - p/E, -2/(alpha + 1)) * L);\
	else var[target] = p*(1 - p)/(pow(1 - p/E, 2)*L);\
    }

void distDNA_F81(unsigned char *x, int n, int s, double *d, double *BF,
		 int variance, double *var, int gamma, double alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p, E;

    L = s;
    E = 1 - BF[0]*BF[0] - BF[1]*BF[1] - BF[2]*BF[2] - BF[3]*BF[3];

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
  	    Nd = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1+= n, s2 += n)
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    COMPUTE_DIST_F81
	    target++;
	}
    }
}

void distDNA_F81_pairdel(unsigned char *x, int n, int s, double *d, double *BF,
			 int variance, double *var, int gamma, double alpha)
{
    int i1, i2, s1, s2, target, Nd, L;
    double p, E;

    E = 1 - BF[0]*BF[0] - BF[1]*BF[1] - BF[2]*BF[2] - BF[3]*BF[3];

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
  	    Nd = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        CHECK_PAIRWISE_DELETION
		if (DifferentBase(x[s1], x[s2])) Nd++;
	    }
	    COMPUTE_DIST_F81
	    target++;
	}
    }
}

#define COUNT_TS_TV1_TV2\
    if (SameBase(x[s1], x[s2])) continue;\
    Nd++;\
    if ((x[s1] | x[s2]) == 152 || (x[s1] | x[s2]) == 104) {\
        Nv1++;\
        continue;\
    }\
    if ((x[s1] | x[s2]) == 168 || (x[s1] | x[s2]) == 88) Nv2++;


#define COMPUTE_DIST_K81\
    P = ((double) (Nd - Nv1 - Nv2)/L);\
    Q = ((double) Nv1/L);\
    R = ((double) Nv2/L);\
    a1 = 1 - 2*P - 2*Q;\
    a2 = 1 - 2*P - 2*R;\
    a3 = 1 - 2*Q - 2*R;\
    d[target] = -0.25*log(a1*a2*a3);\
    if (variance) {\
        a = (1/a1 + 1/a2)/2;\
    	b = (1/a1 + 1/a3)/2;\
    	c = (1/a2 + 1/a3)/2;\
      var[target] = (a*a*P + b*b*Q + c*c*R - pow(a*P + b*Q + c*R, 2))/2;\
    }

void distDNA_K81(unsigned char *x, int n, int s, double *d,
		 int variance, double *var)
{
    int i1, i2, Nd, Nv1, Nv2, L, s1, s2, target;
    double P, Q, R, a1, a2, a3, a, b, c;

    L = s;

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
  	    Nd = Nv1 = Nv2 = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        COUNT_TS_TV1_TV2
	    }
	    COMPUTE_DIST_K81
	    target++;
	}
    }
}

void distDNA_K81_pairdel(unsigned char *x, int n, int s, double *d,
			 int variance, double *var)
{
    int i1, i2, Nd, Nv1, Nv2, L, s1, s2, target;
    double P, Q, R, a1, a2, a3, a, b, c;

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
  	    Nd = Nv1 = Nv2 = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        CHECK_PAIRWISE_DELETION
	        COUNT_TS_TV1_TV2
	    }
	    COMPUTE_DIST_K81
	    target++;
	}
    }
}

#define PREPARE_BF_F84\
    A = (BF[0]*BF[2])/(BF[0] + BF[2]) + (BF[1]*BF[3])/(BF[1] + BF[3]);\
    B = BF[0]*BF[2] + BF[1]*BF[3];\
    C = (BF[0] + BF[2])*(BF[1] + BF[3]);

#define COMPUTE_DIST_F84\
   P = ((double) Ns/L);\
   Q = ((double) (Nd - Ns)/L);\
   d[target] = -2*A*log(1 - P/(2*A) - (A - B)*Q/(2*A*C)) + 2*(A - B - C)*log(1 - Q/(2*C));\
   if (variance) {\
       t1 = A*C;\
       t2 = C*P/2;\
       t3 = (A - B)*Q/2;\
       a = t1/(t1 - t2 - t3);\
       b = A*(A - B)/(t1 - t2 - t3) - (A - B - C)/(C - Q/2);\
       var[target] = (a*a*P + b*b*Q - pow(a*P + b*Q, 2))/L;\
   }

void distDNA_F84(unsigned char *x, int n, int s, double *d,
		 double *BF, int variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, A, B, C, a, b, t1, t2, t3;

    PREPARE_BF_F84
    L = s;

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_F84
	    target++;
	}
    }
}

void distDNA_F84_pairdel(unsigned char *x, int n, int s, double *d,
			 double *BF, int variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, A, B, C, a, b, t1, t2, t3;

    PREPARE_BF_F84

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
		CHECK_PAIRWISE_DELETION
		COUNT_TS_TV
	    }
	    COMPUTE_DIST_F84
	    target++;
	}
    }
}

#define COMPUTE_DIST_T92\
    P = ((double) Ns/L);\
    Q = ((double) (Nd - Ns)/L);\
    a1 = 1 - P/wg - Q;\
    a2 = 1 - 2*Q;\
    d[target] = -wg*log(a1) - 0.5*(1 - wg)*log(a2);\
    if (variance) {\
        c1 = 1/a1;\
        c2 = 1/a2;\
        c3 = wg*(c1 - c2) + c2;\
        var[target] = (c1*c1*P + c3*c3*Q - pow(c1*P + c3*Q, 2))/L;\
    }

void distDNA_T92(unsigned char *x, int n, int s, double *d,
		 double *BF, int variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, wg, a1, a2, c1, c2, c3;

    L = s;
    wg = 2 * (BF[1] + BF[2]) * (1 - (BF[1] + BF[2]));

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_T92
	    target++;
	}
    }
}

void distDNA_T92_pairdel(unsigned char *x, int n, int s, double *d,
			 double *BF, int variance, double *var)
{
    int i1, i2, Nd, Ns, L, target, s1, s2;
    double P, Q, wg, a1, a2, c1, c2, c3;

    wg = 2 * (BF[1] + BF[2]) * (1 - (BF[1] + BF[2]));

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        CHECK_PAIRWISE_DELETION
	        COUNT_TS_TV
	    }
	    COMPUTE_DIST_T92
	    target++;
	}
    }
}

/* returns 1 if one of the base is adenine and
   the other one is guanine surely, 0 otherwise */
#define AdenineAndGuanine(a, b) (a | b) == 200

/* returns 1 if one of the base is cytosine and
   the other one is thymine surely, 0 otherwise */
#define CytosineAndThymine(a, b) (a | b) == 56

#define PREPARE_BF_TN93\
    gR = BF[0] + BF[2];\
    gY = BF[1] + BF[3];\
    k1 = 2 * BF[0] * BF[2] / gR;\
    k2 = 2 * BF[1] * BF[3] / gY;\
    k3 = 2 * (gR * gY - BF[0]*BF[2]*gY/gR - BF[1]*BF[3]*gR/gY);

#define COUNT_TS1_TS2_TV\
    if (DifferentBase(x[s1], x[s2])) {\
        Nd++;\
        if (AdenineAndGuanine(x[s1], x[s2])) {\
            Ns1++;\
    	    continue;\
        }\
        if (CytosineAndThymine(x[s1], x[s2])) Ns2++;\
    }

#define COMPUTE_DIST_TN93\
    P1 = ((double) Ns1/L);\
    P2 = ((double) Ns2/L);\
    Q = ((double) (Nd - Ns1 - Ns2)/L);\
    w1 = 1 - P1/k1 - Q/(2*gR);\
    w2 = 1 - P2/k2 - Q/(2*gY);\
    w3 = 1 - Q/(2*gR*gY);\
    if (variance) {\
        gA2 = BF[0]*BF[0];\
	gC2 = BF[1]*BF[1];\
	gG2 = BF[2]*BF[2];\
	gT2 = BF[3]*BF[3];\
	gAgG = BF[0]*BF[2];\
	gCgT = BF[1]*BF[3];\
	gR2 = gR*gR;\
	gY2 = gY*gY;\
    }\
    if (gamma) {\
        b = -1/alpha;\
	k4 = 2*(BF[0]*BF[2] + BF[1]*BF[3] + gR*gY);\
    	d[target] = alpha * (k1*pow(w1, b) + k2*pow(w2, b) + k3*pow(w3, b) - k4);\
	if (variance) {\
	    b = -(1 + 1/alpha);\
	    c1 = pow(w1, b);\
	    c2 = pow(w2, b);\
    	    c3 = gAgG*c1/gR2 + gCgT*c2/gY2 + ((gA2 + gG2)/(2*gR2) + ((gT2 + gC2)/(2*gY2))) * pow(1 - Q/(2*gR*gY), b);\
	    k4 = c1*P1 + c2*P2 + c3*Q;\
	    var[target] = (c1*c1*P1 + c2*c2*P2 + c3*c3*Q - k4*k4)/L;\
	}\
    } else {\
    	d[target] = -k1*log(w1) - k2*log(w2) - k3*log(w3);\
	if (variance) {\
	    c1 = 1/w1;\
	    c2 = 1/w2;\
	    c3 = 2*gA2*gG2/(gR*(2*gAgG*gR - gR2*P1 - gAgG*Q)) + 2*gC2*gT2/(gY*(2*gCgT*gY - gY2*P2 - gCgT*Q)) + (gR2*(gT2 + gC2) + gY2*(gA2 + gG2))/(2*gR2*gY2 - gR*gY*Q);\
	    k4 = c1*P1 + c2*P2 + c3*Q;\
	    var[target] = (c1*c1*P1 + c2*c2*P2 + c3*c3*Q - k4*k4)/L;\
	}\
    }

void distDNA_TN93(unsigned char *x, int n, int s, double *d,
		  double *BF, int variance, double *var,
		  int gamma, double alpha)
{
    int i1, i2, Nd, Ns1, Ns2, L, target, s1, s2;
    double P1, P2, Q, gR, gY, k1, k2, k3, k4, w1, w2, w3, c1, c2, c3, b;
    double gA2, gC2, gG2, gT2, gAgG, gCgT, gR2, gY2;

    L = s;

    PREPARE_BF_TN93

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns1 = Ns2 = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
		COUNT_TS1_TS2_TV
	    }
	    COMPUTE_DIST_TN93
	    target++;
	}
    }
}

void distDNA_TN93_pairdel(unsigned char *x, int n, int s, double *d,
			  double *BF, int variance, double *var,
			  int gamma, double alpha)
{
    int i1, i2, Nd, Ns1, Ns2, L, target, s1, s2;
    double P1, P2, Q, gR, gY, k1, k2, k3, k4, w1, w2, w3, c1, c2, c3, b;
    double gA2, gC2, gG2, gT2, gAgG, gCgT, gR2, gY2;

    PREPARE_BF_TN93

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns1 = Ns2 = L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
		CHECK_PAIRWISE_DELETION
		COUNT_TS1_TS2_TV
	    }
	    COMPUTE_DIST_TN93
	    target++;
	}
    }
}

void distDNA_GG95(unsigned char *x, int n, int s, double *d,
		  int variance, double *var)
{
    int i1, i2, s1, s2, target, GC, Nd, Ns, tl, npair;
    double *theta, gcprop, *P, pp, *Q, qq, *tstvr, svr, A, sum, ma /* mean alpha */, K1, K2;

    theta = &gcprop;
    P = &pp;
    Q = &qq;
    tstvr = &svr;

    npair = n * (n - 1) / 2;

    theta = (double*)R_alloc(n, sizeof(double));
    P = (double*)R_alloc(npair, sizeof(double));
    Q = (double*)R_alloc(npair, sizeof(double));
    tstvr = (double*)R_alloc(npair, sizeof(double));

    /* get the proportion of GC (= theta) in each sequence */
    for (i1 = 1; i1 <= n; i1++) {
        GC = 0;
	for (s1 = i1 - 1; s1 < i1 + n*(s - 1); s1 += n)
	  if (IsCytosine(x[s1]) || IsGuanine(x[s1])) GC += 1;
	theta[i1 - 1] = ((double) GC / s);
    }

    /* get the proportions of transitions and transversions,
       and the estimates of their ratio for each pair */
    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        COUNT_TS_TV
	    }
	    P[target] = ((double) Ns / s);
	    Q[target] = ((double) (Nd - Ns) / s);
	    A = log(1 - 2*Q[target]);
	    tstvr[target] = 2*(log(1 - 2*P[target] - Q[target]) - 0.5*A)/A;
	    target++;
	}
    }

    /* compute the mean alpha (ma) = mean Ts/Tv */
    sum = 0;
    tl = 0;
    for (i1 = 0; i1 < npair; i1++)
    /* some values of tstvr are -Inf if there is no
       transversions observed: we exclude them */
      if (R_FINITE(tstvr[i1])) {
	  sum += tstvr[i1];
	  tl += 1;
      }
    ma = sum/tl;

    /* compute the distance for each pair */
    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    A = 1 - 2*Q[target];
	    K1 = 1 + ma*(theta[i1 - 1]*(1 - theta[i1 - 1]) + theta[i2 - 1]*(1 - theta[i2 - 1]));
	    K2 = ma*pow(theta[i1 - 1] - theta[i2 - 1], 2)/(ma + 1);
	    d[target] = -0.5*K1*log(A) + K2*(1 - pow(A, 0.25*(ma + 1)));
	    if (variance)
	      var[target] = pow(K1 + K2*0.5*(ma + 1)*pow(A, 0.25*(ma + 1)), 2)*Q[target]*(1 - Q[target])/(A*A * s);
	    target++;
	}
    }
}

void distDNA_GG95_pairdel(unsigned char *x, int n, int s, double *d,
			  int variance, double *var)
{
    int i1, i2, s1, s2, target, *L, length, GC, Nd, Ns, tl, npair;
    double *theta, gcprop, *P, pp, *Q, qq, *tstvr, svr, A, sum, ma /* mean alpha */, K1, K2;

    theta = &gcprop;
    L = &length;
    P = &pp;
    Q = &qq;
    tstvr = &svr;

    npair = n * (n - 1) / 2;

    theta = (double*)R_alloc(n, sizeof(double));
    L = (int*)R_alloc(npair, sizeof(int));
    P = (double*)R_alloc(npair, sizeof(double));
    Q = (double*)R_alloc(npair, sizeof(double));
    tstvr = (double*)R_alloc(npair, sizeof(double));

    /* get the proportion of GC (= theta) in each sequence */
    for (i1 = 1; i1 <= n; i1++) {
        tl = GC = 0;
	for (s1 = i1 - 1; s1 < i1 + n*(s - 1); s1 += n) {
	    if (KnownBase(x[s1])) tl++;
	    else continue;
	    if (IsCytosine(x[s1]) || IsGuanine(x[s1])) GC += 1;
	}
	theta[i1 - 1] = ((double) GC / tl);
    }

    /* get the proportions of transitions and transversions,
       and the estimates of their ratio for each pair; we
       also get the sample size for each pair in L */
    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    Nd = Ns = L[target] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        if (KnownBase(x[s1]) && KnownBase(x[s2])) L[target]++;
		else continue;
	        COUNT_TS_TV
	    }
	    P[target] = ((double) Ns/L[target]);
	    Q[target] = ((double) (Nd - Ns)/L[target]);
	    A = log(1 - 2*Q[target]);
	    tstvr[target] = 2*(log(1 - 2*P[target] - Q[target]) - 0.5*A)/A;
	    target++;
	}
    }

    /* compute the mean alpha (ma) = mean Ts/Tv */
    sum = 0;
    tl = 0;
    for (i1 = 0; i1 < npair; i1++)
    /* some values of tstvr are -Inf if there is no
       transversions observed: we exclude them */
      if (R_FINITE(tstvr[i1])) {
	  sum += tstvr[i1];
	  tl += 1;
      }
    ma = sum/tl;

    /* compute the distance for each pair */
    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    A = 1 - 2*Q[target];
	    K1 = 1 + ma*(theta[i1 - 1]*(1 - theta[i1 - 1]) + theta[i2 - 1]*(1 - theta[i2 - 1]));
	    K2 = ma*pow(theta[i1 - 1] - theta[i2 - 1], 2)/(ma + 1);
	    d[target] = -0.5*K1*log(A) + K2*(1 - pow(A, 0.25*(ma + 1)));
	    if (variance)
	      var[target] = pow(K1 + K2*0.5*(ma + 1)*pow(A, 0.25*(ma + 1)), 2)*Q[target]*(1 - Q[target])/(A*A*L[target]);
	    target++;
	}
    }
}

#define DO_CONTINGENCY_NUCLEOTIDES\
    switch (x[s1]) {\
    case 136 : m = 0; break;\
    case 72 : m = 1; break;\
    case 40 : m = 2; break;\
    case 24 : m = 3; break;\
    }\
    switch (x[s2]) {\
    case 72 : m += 4; break;\
    case 40 : m += 8; break;\
    case 24 : m += 12; break;\
    }\
    Ntab[m]++;

#define COMPUTE_DIST_LogDet\
    for (k = 0; k < 16; k++) Ftab[k] = ((double) Ntab[k]/L);\
    d[target] = -log(detFourByFour(Ftab))/4 - LN4;\
    if (variance) {\
        /* For the inversion, we first make U an identity matrix */\
        for (k = 1; k < 15; k++) U[k] = 0;\
    	U[0] = U[5] = U[10] = U[15] = 1;\
    	/* The matrix is not symmetric, so we use 'dgesv'. */\
    	/* This subroutine puts the result in U. */\
    	F77_CALL(dgesv)(&ndim, &ndim, Ftab, &ndim, ipiv, U, &ndim, &info);\
    	var[target] = (U[0]*U[0]*Ftab[0] + U[1]*U[1]*Ftab[4] +\
    		       U[2]*U[2]*Ftab[8] + U[3]*U[3]*Ftab[12] +\
    		       U[4]*U[4]*Ftab[1] + U[5]*U[5]*Ftab[5] +\
    		       U[6]*U[6]*Ftab[9] + U[7]*U[7]*Ftab[13] +\
    		       U[8]*U[8]*Ftab[2] + U[9]*U[9]*Ftab[6] +\
    		       U[10]*U[10]*Ftab[10] + U[11]*U[11]*Ftab[14] +\
    		       U[12]*U[12]*Ftab[3] + U[13]*U[13]*Ftab[7] +\
    		       U[14]*U[14]*Ftab[11] + U[15]*U[15]*Ftab[15] - 16)/(L*16);\
    }

void distDNA_LogDet(unsigned char *x, int n, int s, double *d,
		    int variance, double *var)
{
    int i1, i2, k, m, s1, s2, target, L, Ntab[16], ndim = 4, info, ipiv[16];
    double Ftab[16], U[16];

    L = s;

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
		DO_CONTINGENCY_NUCLEOTIDES
	    }
	    COMPUTE_DIST_LogDet
	    target++;
	}
    }
}

void distDNA_LogDet_pairdel(unsigned char *x, int n, int s, double *d,
			    int variance, double *var)
{
    int i1, i2, k, m, s1, s2, target, L, Ntab[16], ndim = 4, info, ipiv[16];
    double Ftab[16], U[16];

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    L = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
		CHECK_PAIRWISE_DELETION
		DO_CONTINGENCY_NUCLEOTIDES
	    }
	    COMPUTE_DIST_LogDet
	    target++;
	}
    }
}

void distDNA_BH87(unsigned char *x, int n, int s, double *d)
/* For the moment there is no need to check for pairwise deletions
   since DO_CONTINGENCY_NUCLEOTIDES considers only the known nucleotides.
   In effect the pairwise deletion has possibly been done before.
   The sequence length(s) are used only to compute the variances, which is
   currently not available. */
{
    int i1, i2, k, kb, s1, s2, m, Ntab[16], ROWsums[4];
    double P12[16], P21[16];

    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
		DO_CONTINGENCY_NUCLEOTIDES
	    }

            /* get the rowwise sums of Ntab */
            ROWsums[0] = Ntab[0] + Ntab[4] + Ntab[8] + Ntab[12];
            ROWsums[1] = Ntab[1] + Ntab[5] + Ntab[9] + Ntab[13];
            ROWsums[2] = Ntab[2] + Ntab[6] + Ntab[10] + Ntab[14];
            ROWsums[3] = Ntab[3] + Ntab[7] + Ntab[11] + Ntab[15];

            for (k = 0; k < 16; k++)
              P12[k] = ((double) Ntab[k]);

            /* scale each element of P12 by its rowwise sum */
            for (k = 0; k < 4; k++)
              for (kb = 0; kb < 16; kb += 4)
            	P12[k + kb] = P12[k + kb]/ROWsums[k];

            d[n*(i2 - 1) + i1 - 1] = -log(detFourByFour(P12))/4;

            /* compute the columnwise sums of Ntab: these
               are the rowwise sums of its transpose */
            ROWsums[0] = Ntab[0] + Ntab[1] + Ntab[2] + Ntab[3];
            ROWsums[1] = Ntab[4] + Ntab[5] + Ntab[6] + Ntab[7];
            ROWsums[2] = Ntab[8] + Ntab[9] + Ntab[10] + Ntab[11];
            ROWsums[3] = Ntab[12] + Ntab[13] + Ntab[14] + Ntab[15];

            /* transpose Ntab and store the result in P21 */
            for (k = 0; k < 4; k++)
               for (kb = 0; kb < 4; kb++)
            	 P21[kb + 4*k] = Ntab[k + 4*kb];

            /* scale as above */
            for (k = 0; k < 4; k++)
              for (kb = 0; kb < 16; kb += 4)
            	P21[k + kb] = P21[k + kb]/ROWsums[k];

            d[n*(i1 - 1) + i2 - 1] = -log(detFourByFour(P21))/4;
	}
    }
}

#define COMPUTE_DIST_ParaLin\
    for (k = 0; k < 16; k++) Ftab[k] = ((double) Ntab[k]/L);\
    d[target] = -log(detFourByFour(Ftab)/\
		     sqrt(find[0][i1 - 1]*find[1][i1 - 1]*find[2][i1 - 1]*find[3][i1 - 1]*\
			  find[0][i2 - 1]*find[1][i2 - 1]*find[2][i2 - 1]*find[3][i2 - 1]))/4;\
    if (variance) {\
        /* For the inversion, we first make U an identity matrix */\
        for (k = 1; k < 15; k++) U[k] = 0;\
    	U[0] = U[5] = U[10] = U[15] = 1;\
    	/* The matrix is not symmetric, so we use 'dgesv'. */\
    	/* This subroutine puts the result in U. */\
    	F77_CALL(dgesv)(&ndim, &ndim, Ftab, &ndim, ipiv, U, &ndim, &info);\
    	var[target] = (U[0]*U[0]*Ftab[0] + U[1]*U[1]*Ftab[4] +\
    		       U[2]*U[2]*Ftab[8] + U[3]*U[3]*Ftab[12] +\
    		       U[4]*U[4]*Ftab[1] + U[5]*U[5]*Ftab[5] +\
    		       U[6]*U[6]*Ftab[9] + U[7]*U[7]*Ftab[13] +\
    		       U[8]*U[8]*Ftab[2] + U[9]*U[9]*Ftab[6] +\
    		       U[10]*U[10]*Ftab[10] + U[11]*U[11]*Ftab[14] +\
    		       U[12]*U[12]*Ftab[3] + U[13]*U[13]*Ftab[7] +\
    		       U[14]*U[14]*Ftab[11] + U[15]*U[15]*Ftab[15] -\
		       4*(1/sqrt(find[0][i1 - 1]*find[0][i2 - 1]) +\
                       1/sqrt(find[1][i1 - 1]*find[1][i2 - 1]) +\
		       1/sqrt(find[2][i1 - 1]*find[2][i2 - 1]) +\
                       1/sqrt(find[3][i1 - 1]*find[3][i2 - 1])))/(L*16);\
    }

void distDNA_ParaLin(unsigned char *x, int n, int s, double *d,
		     int variance, double *var)
{
    int i1, i2, k, s1, s2, m, target, L, Ntab[16], ndim = 4, info, ipiv[16];
    double Ftab[16], U[16], *find[4];

    L = s;

    for (k = 0; k < 4; k++)
      find[k] = (double*)R_alloc(n, sizeof(double));

    for (i1 = 0; i1 < n; i1++)
      for (k = 0; k < 4; k++) find[k][i1] = 0.0;

    for (i1 = 0; i1 < n; i1++) {
        for (s1 = i1; s1 < i1 + n*(s - 1) + 1; s1+= n) {
            switch (x[s1]) {
	    case 136 : find[0][i1]++; break;
	    case 40 : find[1][i1]++; break;
	    case 72 : find[2][i1]++; break;
	    case 24 : find[3][i1]++; break;
	    }
        }
        for (k = 0; k < 4; k++) find[k][i1] /= L;
    }

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
		DO_CONTINGENCY_NUCLEOTIDES
	    }
	    COMPUTE_DIST_ParaLin
	    target++;
	}
    }
}

void distDNA_ParaLin_pairdel(unsigned char *x, int n, int s, double *d,
			     int variance, double *var)
{
    int i1, i2, k, s1, s2, m, target, L, Ntab[16], ndim = 4, info, ipiv[16];
    double Ftab[16], U[16], *find[4];

    L = 0;

    for (k = 0; k < 4; k++)
      find[k] = (double*)R_alloc(n, sizeof(double));

    for (i1 = 0; i1 < n; i1++)
      for (k = 0; k < 4; k++) find[k][i1] = 0.0;

    for (i1 = 0; i1 < n; i1++) {
        L = 0;
        for (s1 = i1; s1 < i1 + n*(s - 1) + 1; s1+= n) {
	    if (KnownBase(x[s1])) {
	        L++;
                switch (x[s1]) {
	        case 136 : find[0][i1]++; break;
	        case 40 : find[1][i1]++; break;
	        case 72 : find[2][i1]++; break;
	        case 24 : find[3][i1]++; break;
	        }
	    }
        }
        for (k = 0; k < 4; k++) find[k][i1] /= L;
    }

    target = 0;
    for (i1 = 1; i1 < n; i1++) {
        for (i2 = i1 + 1; i2 <= n; i2++) {
	    L = 0;
	    for (k = 0; k < 16; k++) Ntab[k] = 0;
	    for (s1 = i1 - 1, s2 = i2 - 1; s1 < i1 + n*(s - 1); s1 += n, s2 += n) {
	        CHECK_PAIRWISE_DELETION
		DO_CONTINGENCY_NUCLEOTIDES
	    }
	    COMPUTE_DIST_ParaLin
	    target++;
	}
    }
}

/* a look-up table is much faster than switch (2012-01-10) */
SEXP BaseProportion(SEXP x)
{
    long i;
    unsigned char *p;
    double n, count[256], *BF;
    SEXP res;

    PROTECT(x = coerceVector(x, RAWSXP));
    memset(count, 0, 256*sizeof(double));
    n = XLENGTH(x);
    p = RAW(x);

    for (i = 0; i < n; i++) count[p[i]]++;

    PROTECT(res = allocVector(REALSXP, 17));
    BF = REAL(res);
    BF[0] = count[136];
    BF[1] = count[40];
    BF[2] = count[72];
    BF[3] = count[24];
    BF[4] = count[192];
    BF[5] = count[160];
    BF[6] = count[144];
    BF[7] = count[96];
    BF[8] = count[80];
    BF[9] = count[48];
    BF[10] = count[224];
    BF[11] = count[176];
    BF[12] = count[208];
    BF[13] = count[112];
    BF[14] = count[240];
    BF[15] = count[4];
    BF[16] = count[2];
    UNPROTECT(2);
    return res;
}

#define SEGCOL seg[j] = 1; done = 1; break

void seg_sites_a(unsigned char *x, int *seg, int n, int s)
{
    long i, end, done;
    int j;
    unsigned char base;

    for (j = 0; j < s; j++) {

        i = (long) n * j; /* start */
	end = i + n - 1;

        base = x[i];
	done = 0;

	while (!KnownBase(base)) {
	    /* in this while-loop, we are not yet sure that 'base' is known,
	       so we must be careful with the comparisons */
	    i++;
	    if (i > end) {
		done = 1;
		break;
	    }
	    if (base != x[i]) {
		if (base != 2 && x[i] != 2) { /* both should not be "?" */
		    if (base > 4) {
			if (x[i] == 4) { /* 'base' is not a gap but x[i] is one => this is a segregating site */
			    SEGCOL;
			} else { /* both are an ambiguous base */
			    if (DifferentBase(x[i], base)) {
				SEGCOL;
			    }
			}
		    } else { /* 'base' is a gap but x[i] is different => this is a segregating site */
			SEGCOL;
		    }
		}
		base = x[i];
	    }
	}

	if (done) continue;

	i++;
	while (i <= end) {
	    if (x[i] != base) {
		if (x[i] == 4) {
		    SEGCOL;
		} else {
		    if (DifferentBase(x[i], base)) {
			SEGCOL;
		    }
		}
	    }
	    i++;
	}
    }
}

void seg_sites_strict(unsigned char *x, int *seg, int n, int s)
{
    long i, end;
    int j;
    unsigned char b;

    for (j = 0; j < s; j++) {
	i = (long) n * j; /* start */
	end = i + n - 1;
        b = x[i];
	i++;
	while (i <= end) {
	    if (x[i] != b) {
		seg[j] = 1;
		break;
	    }
	    i++;
	}
    }
}

SEXP SegSites(SEXP DNASEQ, SEXP STRICT)
{
    int n, s, *seg;
    unsigned char *x;
    SEXP ans;

    PROTECT(STRICT = coerceVector(STRICT, INTSXP));
    PROTECT(DNASEQ = coerceVector(DNASEQ, RAWSXP));
    x = RAW(DNASEQ);
    n = nrows(DNASEQ);
    s = ncols(DNASEQ);

    PROTECT(ans = allocVector(INTSXP, s));
    seg = INTEGER(ans);
    memset(seg, 0, s * sizeof(int));

    if (INTEGER(STRICT)[0]) {
	seg_sites_strict(x, seg, n, s);
    } else {
	seg_sites_a(x, seg, n, s);
    }

    UNPROTECT(3);
    return ans;
}

SEXP GlobalDeletionDNA(SEXP DNASEQ)
{
    int i, j, n, s;
    unsigned char *x;
    int *keep;
    SEXP res;

    PROTECT(DNASEQ = coerceVector(DNASEQ, RAWSXP));
    x = RAW(DNASEQ);
    n = nrows(DNASEQ);
    s = ncols(DNASEQ);

    PROTECT(res = allocVector(INTSXP, s));
    keep = INTEGER(res);
    memset(keep, 1, s * sizeof(int));

    for (j = 0; j < s; j++) {
        i = n * j;
	while (i < n * (j + 1)) {
	    if (KnownBase(x[i])) i++;
	    else {
	        keep[j] = 0;
		break;
	    }
	}
    }
    UNPROTECT(2);
    return res;
}

SEXP dist_dna(SEXP DNASEQ, SEXP MODEL, SEXP BASEFREQ, SEXP PAIRDEL,
	      SEXP VARIANCE, SEXP GAMMA, SEXP ALPHA)
{
    int n, s, model, pairdel, variance, gamma, Ndist;
    double *BF, alpha, *d, *var;
    unsigned char *x;
    SEXP res, distvar;

    PROTECT(DNASEQ = coerceVector(DNASEQ, RAWSXP));
    PROTECT(MODEL = coerceVector(MODEL, INTSXP));
    PROTECT(BASEFREQ = coerceVector(BASEFREQ, REALSXP));
    PROTECT(PAIRDEL = coerceVector(PAIRDEL, INTSXP));
    PROTECT(VARIANCE = coerceVector(VARIANCE, INTSXP));
    PROTECT(GAMMA = coerceVector(GAMMA, INTSXP));
    PROTECT(ALPHA = coerceVector(ALPHA, REALSXP));

    n = nrows(DNASEQ);
    s = ncols(DNASEQ);
    x = RAW(DNASEQ);
    model = INTEGER(MODEL)[0];
    Ndist = n * (n - 1) / 2;
    if (model == 11) Ndist = n * n;
    BF = REAL(BASEFREQ);
    pairdel = INTEGER(PAIRDEL)[0];
    variance = INTEGER(VARIANCE)[0];
    if (variance) {
	PROTECT(distvar =  allocVector(REALSXP, Ndist));
	var = REAL(distvar);
    }
    gamma = INTEGER(GAMMA)[0];
    if (gamma) alpha = REAL(ALPHA)[0];

    PROTECT(res =  allocVector(REALSXP, Ndist));
    d = REAL(res);

    switch (model) {
    case 1 :
	if (pairdel) {
	    distDNA_raw_pairdel(x, n, s, d, 1);
	} else {
	    distDNA_raw(x, n, s, d, 1);
	}
	break;
    case 2 :
	if (pairdel) {
	    distDNA_JC69_pairdel(x, n, s, d, variance, var, gamma, alpha);
	} else {
	    distDNA_JC69(x, n, s, d, variance, var, gamma, alpha);
	}
	break;
    case 3 :
	if (pairdel) {
	    distDNA_K80_pairdel(x, n, s, d, variance, var, gamma, alpha);
	} else {
	    distDNA_K80(x, n, s, d, variance, var, gamma, alpha);
	}
	break;
    case 4 :
	if (pairdel) {
	    distDNA_F81_pairdel(x, n, s, d, BF, variance, var, gamma, alpha);
	} else {
	    distDNA_F81(x, n, s, d, BF, variance, var, gamma, alpha);
	}
	break;
    case 5 :
	if (pairdel) {
	    distDNA_K81_pairdel(x, n, s, d, variance, var);
	} else {
	    distDNA_K81(x, n, s, d, variance, var);
	}
	break;
    case 6 :
	if (pairdel) {
	    distDNA_F84_pairdel(x, n, s, d, BF, variance, var);
	} else {
	    distDNA_F84(x, n, s, d, BF, variance, var);
	}
	break;
    case 7 :
	if (pairdel) {
	    distDNA_T92_pairdel(x, n, s, d, BF, variance, var);
	} else {
	    distDNA_T92(x, n, s, d, BF, variance, var);
	}
	break;
    case 8 :
	if (pairdel) {
	    distDNA_TN93_pairdel(x, n, s, d, BF, variance, var, gamma, alpha);
	} else {
	    distDNA_TN93(x, n, s, d, BF, variance, var, gamma, alpha);
	}
	break;
    case 9 :
	if (pairdel) {
	    distDNA_GG95_pairdel(x, n, s, d, variance, var);
	} else {
	    distDNA_GG95(x, n, s, d, variance, var);
	}
	break;
    case 10 :
	if (pairdel) {
	    distDNA_LogDet_pairdel(x, n, s, d, variance, var);
	} else {
	    distDNA_LogDet(x, n, s, d, variance, var);
	}
	break;
    case 11 :
	distDNA_BH87(x, n, s, d);
	break;
    case 12 :
	if (pairdel) {
	    distDNA_ParaLin_pairdel(x, n, s, d, variance, var);
	} else {
	    distDNA_ParaLin(x, n, s, d, variance, var);
	}
	break;
    case 13 :
	if (pairdel) {
	    distDNA_raw_pairdel(x, n, s, d, 0);
	} else {
	    distDNA_raw(x, n, s, d, 0);
	}
	break;
    case 14 :
	if (pairdel) {
	    distDNA_TsTv(x, n, s, d, 1, 1);
	} else {
	    distDNA_TsTv(x, n, s, d, 1, 0);
	}
	break;
    case 15 :
	if (pairdel) {
	    distDNA_TsTv(x, n, s, d, 0, 1);
	} else {
	    distDNA_TsTv(x, n, s, d, 0, 0);
	}
	break;
    case 16 :
	distDNA_indel(x, n, s, d);
	break;
    case 17 :
	distDNA_indelblock(x, n, s, d);
	break;
    }

    if (variance) {
	SEXP obj;
	PROTECT(obj = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(obj, 0, res);
	SET_VECTOR_ELT(obj, 1, distvar);
	UNPROTECT(10);
	return obj;
    } else {
	UNPROTECT(8);
	return res;
    }
}

SEXP C_where(SEXP DNASEQ, SEXP PAT)
{
	int p, j, nans;
	double s, *buf, *a;
	long i, k;
	unsigned char *x, *pat;
	SEXP ans;

	PROTECT(DNASEQ = coerceVector(DNASEQ, RAWSXP));
	PROTECT(PAT = coerceVector(PAT, RAWSXP));
	x = RAW(DNASEQ);
	pat = RAW(PAT);
	s = XLENGTH(DNASEQ);
	p = LENGTH(PAT);
	nans = 0;

	buf = (double *)R_alloc(s, sizeof(double));

	for (i = 0; i <= s - p; i++) {
		k = i; j = 0;
		while (1) {
			if (x[k] != pat[j]) break;
			j++; k++;
			if (j == p) {
				buf[nans++] = i + 1;
				break;
			}
		}
	}

	PROTECT(ans = allocVector(REALSXP, nans));
	if (nans) {
	    a = REAL(ans);
	    for (i = 0; i < nans; i++) a[i] = buf[i];
	}

	UNPROTECT(3);
	return ans;
}

unsigned char codon2aa_Code1(unsigned char x, unsigned char y, unsigned char z)
{
    if (KnownBase(x)) {
	if (IsAdenine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x4b; /* codon is AAR => 'K' */
		    if (IsPyrimidine(z)) return 0x4e; /* codon is AAY => 'N' */
		    return 0x58; /* 'X' */
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x54; /* codon is ACN => 'T' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsPurine(z)) return 0x52; /* codon is AGR => 'R' */
		    if (IsPyrimidine(z)) return 0x53; /* codon is AGY => 'S' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsGuanine(z)) return 0x4d; /* codon is ATG => 'M' */
		    if (z & 176) return 0x49; /* codon is ATH => 'I' */
		    return 0x58;
		}
	    }
	    return 0x58;
	}
	if (IsCytosine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x51; /* codon is CAR => 'Q'*/
		if (IsPyrimidine(z)) return 0x48; /* codon is CAY => 'H' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x50; /* codon is CCN => 'P'*/
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x52; /* codon is CGN => 'R' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x4c; /* codon is CTN => 'L' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsGuanine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x45; /* codon is GAR => 'E' */
		if (IsPyrimidine(z)) return 0x44; /* codon is GAY => 'D' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x41; /* codon is GCN => 'A' */
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x47; /* codon is GGN => 'G' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x56; /* codon is GTN => 'V' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsThymine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x2a; /* codon is TAR => '*' */
		    if (IsPyrimidine(z)) return 0x59; /* codon is TAY => 'Y' */
		    return 0x58;
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x53; /* codon is TCN => 'S' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsAdenine(z)) return 0x2a; /* codon is TGA => '*' */
		    if (IsGuanine(z)) return 0x57; /* codon is TGG => 'W' */
		    if (IsPyrimidine(z)) return 0x43; /* codon is TGY => 'C' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsPurine(z)) return 0x4c; /* codon is TTR => 'L' */
		    if (IsPyrimidine(z)) return 0x46; /* codon is TTY => 'F' */
		    return 0x58;
		}
	    } else if (IsPurine(y) & IsAdenine(z)) return 0x2a; /* codon is TRA => '*' */
	    return 0x58;
	}
    } else {
	if ((x == 144) && IsThymine(y) && IsPurine(z)) return 0x52; /* codon is MGR => 'R'*/
	if ((x == 48) && IsThymine(y) && IsPurine(z)) return 0x4c; /* codon is YTR => 'L'*/
    }
    return 0x58;
}

unsigned char codon2aa_Code2(unsigned char x, unsigned char y, unsigned char z)
{
    if (KnownBase(x)) {
	if (IsAdenine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x4b; /* codon is AAR => 'K' */
		    if (IsPyrimidine(z)) return 0x4e; /* codon is AAY => 'N' */
		    return 0x58; /* 'X' */
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x54; /* codon is ACN => 'T' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsPurine(z)) return 0x2a; /* codon is AGR => '*' */
		    if (IsPyrimidine(z)) return 0x53; /* codon is AGY => 'S' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsPurine(z)) return 0x4d; /* codon is ATR => 'M' */
		    if (IsPyrimidine(z)) return 0x49; /* codon is ATY => 'I' */
		    return 0x58;
		}
	    }
	    return 0x58;
	}
	if (IsCytosine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x51; /* codon is CAR => 'Q'*/
		if (IsPyrimidine(z)) return 0x48; /* codon is CAY => 'H' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x50; /* codon is CCN => 'P'*/
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x52; /* codon is CGN => 'R' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x4c; /* codon is CTN => 'L' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsGuanine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x45; /* codon is GAR => 'E' */
		if (IsPyrimidine(z)) return 0x44; /* codon is GAY => 'D' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x41; /* codon is GCN => 'A' */
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x47; /* codon is GGN => 'G' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x56; /* codon is GTN => 'V' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsThymine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x2a; /* codon is TAR => '*' */
		    if (IsPyrimidine(z)) return 0x59; /* codon is TAY => 'Y' */
		    return 0x58;
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x53; /* codon is TCN => 'S' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsPurine(z)) return 0x57; /* codon is TGR => 'W' */
		    if (IsPyrimidine(z)) return 0x43; /* codon is TGY => 'C' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsPurine(z)) return 0x4c; /* codon is TTR => 'L' */
		    if (IsPyrimidine(z)) return 0x46; /* codon is TTY => 'F' */
		    return 0x58;
		}
	    }
	    return 0x58;
	}
    } else {
	if ((x == 48) && IsThymine(y) && IsPurine(z)) return 0x4c; /* codon is YTR => 'L'*/
    }
    return 0x58;
}

unsigned char codon2aa_Code3(unsigned char x, unsigned char y, unsigned char z)
{
    if (KnownBase(x)) {
	if (IsAdenine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x4b; /* codon is AAR => 'K' */
		    if (IsPyrimidine(z)) return 0x4e; /* codon is AAY => 'N' */
		    return 0x58; /* 'X' */
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x54; /* codon is ACN => 'T' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsPurine(z)) return 0x52; /* codon is AGR => 'R' */
		    if (IsPyrimidine(z)) return 0x53; /* codon is AGY => 'S' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsPurine(z)) return 0x4d; /* codon is ATR => 'M' */
		    if (IsPyrimidine(z)) return 0x49; /* codon is ATY => 'I' */
		    return 0x58;
		}
	    }
	    return 0x58;
	}
	if (IsCytosine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x51; /* codon is CAR => 'Q'*/
		if (IsPyrimidine(z)) return 0x48; /* codon is CAY => 'H' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x50; /* codon is CCN => 'P'*/
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x52; /* codon is CGN => 'R' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x4c; /* codon is CTN => 'T' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsGuanine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x45; /* codon is GAR => 'E' */
		if (IsPyrimidine(z)) return 0x44; /* codon is GAY => 'D' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x41; /* codon is GCN => 'A' */
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x47; /* codon is GGN => 'G' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x56; /* codon is GTN => 'V' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsThymine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x2a; /* codon is TAR => '*' */
		    if (IsPyrimidine(z)) return 0x59; /* codon is TAY => 'Y' */
		    return 0x58;
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x53; /* codon is TCN => 'S' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsPurine(z)) return 0x57; /* codon is TGR => 'W' */
		    if (IsPyrimidine(z)) return 0x43; /* codon is TGY => 'C' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsPurine(z)) return 0x4c; /* codon is TTR => 'L' */
		    if (IsPyrimidine(z)) return 0x46; /* codon is TTY => 'F' */
		    return 0x58;
		}
	    } else if (IsPurine(y) & IsAdenine(z)) return 0x2a; /* codon is TRA => '*' */
	    return 0x58;
	}
    } else {
	if ((x == 144) && IsThymine(y) && IsPurine(z)) return 0x52; /* codon is MGR => 'R'*/
	if ((x == 48) && IsThymine(y) && IsPurine(z)) return 0x4c; /* codon is YTR => 'L'*/
    }
    return 0x58;
}

unsigned char codon2aa_Code4(unsigned char x, unsigned char y, unsigned char z)
{
    if (KnownBase(x)) {
	if (IsAdenine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x4b; /* codon is AAR => 'K' */
		    if (IsPyrimidine(z)) return 0x4e; /* codon is AAY => 'N' */
		    return 0x58; /* 'X' */
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x54; /* codon is ACN => 'T' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsPurine(z)) return 0x52; /* codon is AGR => 'R' */
		    if (IsPyrimidine(z)) return 0x53; /* codon is AGY => 'S' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsGuanine(z)) return 0x4d; /* codon is ATG => 'M' */
		    if (z & 176) return 0x49; /* codon is ATH => 'I' */
		    return 0x58;
		}
	    }
	    return 0x58;
	}
	if (IsCytosine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x51; /* codon is CAR => 'Q'*/
		if (IsPyrimidine(z)) return 0x48; /* codon is CAY => 'H' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x50; /* codon is CCN => 'P'*/
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x52; /* codon is CGN => 'R' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x4c; /* codon is CTN => 'L' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsGuanine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x45; /* codon is GAR => 'E' */
		if (IsPyrimidine(z)) return 0x44; /* codon is GAY => 'D' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x41; /* codon is GCN => 'A' */
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x47; /* codon is GGN => 'G' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x56; /* codon is GTN => 'V' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsThymine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x2a; /* codon is TAR => '*' */
		    if (IsPyrimidine(z)) return 0x59; /* codon is TAY => 'Y' */
		    return 0x58;
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x53; /* codon is TCN => 'S' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsPurine(z)) return 0x57; /* codon is TGR => 'W' */
		    if (IsPyrimidine(z)) return 0x43; /* codon is TGY => 'C' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsPurine(z)) return 0x4c; /* codon is TTR => 'L' */
		    if (IsPyrimidine(z)) return 0x46; /* codon is TTY => 'F' */
		    return 0x58;
		}
	    } else if (IsPurine(y) & IsAdenine(z)) return 0x2a; /* codon is TRA => '*' */
	    return 0x58;
	}
    } else {
	if ((x == 144) && IsThymine(y) && IsPurine(z)) return 0x52; /* codon is MGR => 'R'*/
	if ((x == 48) && IsThymine(y) && IsPurine(z)) return 0x4c; /* codon is YTR => 'L'*/
    }
    return 0x58;
}

unsigned char codon2aa_Code5(unsigned char x, unsigned char y, unsigned char z)
{
    if (KnownBase(x)) {
	if (IsAdenine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x4b; /* codon is AAR => 'K' */
		    if (IsPyrimidine(z)) return 0x4e; /* codon is AAY => 'N' */
		    return 0x58; /* 'X' */
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x54; /* codon is ACN => 'T' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (z > 4) return 0x53; /* codon is AGN => 'S' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsPurine(z)) return 0x4d; /* codon is ATR => 'M' */
		    if (IsPyrimidine(z)) return 0x49; /* codon is ATY => 'I' */
		    return 0x58;
		}
	    }
	    return 0x58;
	}
	if (IsCytosine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x51; /* codon is CAR => 'Q'*/
		if (IsPyrimidine(z)) return 0x48; /* codon is CAY => 'H' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x50; /* codon is CCN => 'P'*/
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x52; /* codon is CGN => 'R' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x4c; /* codon is CTN => 'L' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsGuanine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x45; /* codon is GAR => 'E' */
		if (IsPyrimidine(z)) return 0x44; /* codon is GAY => 'D' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x41; /* codon is GCN => 'A' */
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x47; /* codon is GGN => 'G' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x56; /* codon is GTN => 'V' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsThymine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x2a; /* codon is TAR => '*' */
		    if (IsPyrimidine(z)) return 0x59; /* codon is TAY => 'Y' */
		    return 0x58;
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x53; /* codon is TCN => 'S' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsPurine(z)) return 0x57; /* codon is TGR => 'W' */
		    if (IsPyrimidine(z)) return 0x43; /* codon is TGY => 'C' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsPurine(z)) return 0x4c; /* codon is TTR => 'L' */
		    if (IsPyrimidine(z)) return 0x46; /* codon is TTY => 'F' */
		    return 0x58;
		}
	    } else if (IsPurine(y) & IsAdenine(z)) return 0x2a; /* codon is TRA => '*' */
	    return 0x58;
	}
    } else {
	if ((x == 144) && IsThymine(y) && IsPurine(z)) return 0x52; /* codon is MGR => 'R'*/
	if ((x == 48) && IsThymine(y) && IsPurine(z)) return 0x4c; /* codon is YTR => 'L'*/
    }
    return 0x58;
}

unsigned char codon2aa_Code6(unsigned char x, unsigned char y, unsigned char z)
{
    if (KnownBase(x)) {
	if (IsAdenine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x4b; /* codon is AAR => 'K' */
		    if (IsPyrimidine(z)) return 0x4e; /* codon is AAY => 'N' */
		    return 0x58; /* 'X' */
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x54; /* codon is ACN => 'T' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsPurine(z)) return 0x52; /* codon is AGR => 'R' */
		    if (IsPyrimidine(z)) return 0x53; /* codon is AGY => 'S' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsGuanine(z)) return 0x4d; /* codon is ATG => 'M' */
		    if (z & 176) return 0x49; /* codon is ATH => 'I' */
		    return 0x58;
		}
	    }
	    return 0x58;
	}
	if (IsCytosine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x51; /* codon is CAR => 'Q'*/
		if (IsPyrimidine(z)) return 0x48; /* codon is CAY => 'H' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x50; /* codon is CCN => 'P'*/
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x52; /* codon is CGN => 'R' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x4c; /* codon is CTN => 'L' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsGuanine(x)) {
	    if (IsAdenine(y)) {
		if (IsPurine(z)) return 0x45; /* codon is GAR => 'E' */
		if (IsPyrimidine(z)) return 0x44; /* codon is GAY => 'D' */
		return 0x58;
	    }
	    if (IsCytosine(y)) {
		if (z > 4) return 0x41; /* codon is GCN => 'A' */
		return 0x58;
	    }
	    if (IsGuanine(y)) {
		if (z > 4) return 0x47; /* codon is GGN => 'G' */
		return 0x58;
	    }
	    if (IsThymine(y)) {
		if (z > 4) return 0x56; /* codon is GTN => 'V' */
		return 0x58;
	    }
	    return 0x58;
	}
	if (IsThymine(x)) {
	    if (KnownBase(y)) {
		if (IsAdenine(y)) {
		    if (IsPurine(z)) return 0x2a; /* codon is TAR => 'Q' */
		    if (IsPyrimidine(z)) return 0x59; /* codon is TAY => 'Y' */
		    return 0x58;
		}
		if (IsCytosine(y)) {
		    if (z > 4) return 0x53; /* codon is TCN => 'S' */
		    return 0x58;
		}
		if (IsGuanine(y)) {
		    if (IsAdenine(z)) return 0x2a; /* codon is TGA => '*' */
		    if (IsGuanine(z)) return 0x57; /* codon is TGG => 'W' */
		    if (IsPyrimidine(z)) return 0x43; /* codon is TGY => 'C' */
		    return 0x58;
		}
		if (IsThymine(y)) {
		    if (IsPurine(z)) return 0x4c; /* codon is TTR => 'L' */
		    if (IsPyrimidine(z)) return 0x46; /* codon is TTY => 'F' */
		    return 0x58;
		}
	    } else if (IsPurine(y) & IsAdenine(z)) return 0x2a; /* codon is TRA => '*' */
	    return 0x58;
	}
    } else {
	if ((x == 144) && IsThymine(y) && IsPurine(z)) return 0x52; /* codon is MGR => 'R'*/
	if ((x == 48) && IsThymine(y) && IsPurine(z)) return 0x4c; /* codon is YTR => 'L'*/
    }
    return 0x58;
}

void trans_DNA2AA(unsigned char *x, int *s, unsigned char *res, int *code)
{
    int i = 0, j = 0;
    unsigned char (*FUN)(unsigned char x, unsigned char y, unsigned char z);

    /* NOTE: using 'switch' provokes a memory leak */
    while (1) {
	if (*code == 1) { FUN = &codon2aa_Code1; break; }
	if (*code == 2) { FUN = &codon2aa_Code2; break; }
	if (*code == 3) { FUN = &codon2aa_Code3; break; }
	if (*code == 4) { FUN = &codon2aa_Code4; break; }
	if (*code == 5) { FUN = &codon2aa_Code5; break; }
	if (*code == 6) { FUN = &codon2aa_Code6; break; }
    }

    while (i < *s) {
	res[j] = FUN(x[i], x[i + 1], x[i + 2]);
	j++; i += 3;
    }
}

SEXP leading_trailing_gaps_to_N(SEXP DNASEQ)
{
    int i, n, s;
    long j, k;
    unsigned char *x, *z;
    SEXP ans;

    PROTECT(DNASEQ = coerceVector(DNASEQ, RAWSXP));
    x = RAW(DNASEQ);
    n = nrows(DNASEQ);
    s = ncols(DNASEQ);

    PROTECT(ans = allocVector(RAWSXP, n * s));
    z = RAW(ans);

    memcpy(z, x, n * s);

    for (i = 0; i < n; i++) { /* leading gaps */
	j = (long) i; /* start of the seq */
	k = (long) j + n * (s - 1); /* last site of seq j */
	while (x[j] == 4 && j <= k) {
	    z[j] = 240; /* - -> N */
	    j += n;
	}
    }
    for (i = 0; i < n; i++) { /* trailing gaps */
	k = (long) i; /* start of the seq */
	j = (long) k + n * (s - 1); /* the last site of seq j */
	while (x[j] == 4 && j >= k) {
		z[j] = 240; /* - -> N */
		j -= n;
	}
    }

    UNPROTECT(2);
    return ans;
}
