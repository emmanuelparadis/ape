/* mvr.c    2012-05-02 */

/* Copyright 2011-2012 Andrei-Alin Popescu */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

void C_mvr(double *D, double* v,int *N, int *edge1, int *edge2, double *edge_length)
{
	double *S, Sdist, *new_v, Ndist, *new_dist, A, B, smallest_S;
	int n, i, j, k, ij, smallest, OTU1, OTU2, cur_nod, o_l, *otu_label;

	S = &Sdist;
	new_dist = &Ndist;
	otu_label = &o_l;
        n = *N;
	cur_nod = 2*n - 2;

	S = (double*)R_alloc(n + 1, sizeof(double));
	new_dist = (double*)R_alloc(n*(n - 1)/2, sizeof(double));
        new_v = (double*)R_alloc(n*(n - 1)/2, sizeof(double));
	otu_label = (int*)R_alloc(n + 1, sizeof(int));

	for (i = 1; i <= n; i++) otu_label[i] = i; /* otu_label[0] is not used */

	k = 0;

	while (n > 3) {

                /*for(i=1;i<n;i++)
                  {
                    for(j=i+1;j<=n;j++)
                      {
                        Rprintf("d[%i,%i]=%f ",i,j,D[give_index(i,j,n)]);
                      }
                    Rprintf("\n");
                  } */

		for (i = 1; i <= n; i++)
                {double sum=0;
                  for( j=1;j<=n;j++)
                   {if(i==j)continue;
                     //Rprintf("index(%i,%i)=%i\n",i,j,give_index(i,j,n));
                     //Rprintf("D[%i,%i]=%f\n",i,j,D[give_index(i,j,n)]);
                     sum+=D[give_index(i,j,n)];
                   }
                S[i]=sum;
                //Rprintf("\n");
                //Rprintf("S[%i]=%f\n",i,S[i]);
                //Rprintf("\n");
                }
		ij = 0;
		smallest_S = 1e50;
		B = n - 2;
		for (i = 1; i < n; i++) {
			for (j = i + 1; j <= n; j++) {

				A = B*D[ij] - S[i] - S[j];
                                /*Rprintf("D[ij]=%f\n",D[ij]);
                                Rprintf("S[%i]=%f\n",i,S[i]);
                                Rprintf("S[%i]=%f\n",j,S[j]);
                                Rprintf("B=%f\n",B);
                                Rprintf("A=%f\n",(B*D[ij] - S[i] - S[j]));
                                Rprintf("Q[%i,%i]=%f\n",i,j,A);*/
				if (A < smallest_S) {
					OTU1 = i;
					OTU2 = j;
					smallest_S = A;
					smallest = ij;
				}
				ij++;
			}
		}

                //Rprintf("agglomerating %i and %i, Q=%f \n",OTU1,OTU2,smallest_S);

                /*for(i=1;i<n;i++)
                  {
                    for(j=i+1;j<=n;j++)
                      {
                        Rprintf("d[%i,%i]=%f ",i,j,D[give_index(i,j,n)]);
                      }
                    Rprintf("\n");
                  }

                for(i=1;i<n;i++)
                  {
                    for(j=i+1;j<=n;j++)
                      {
                        Rprintf("v[%i,%i]=%f ",i,j,v[give_index(i,j,n)]);
                      }
                    Rprintf("\n");
                  }*/

		edge2[k] = otu_label[OTU1];
		edge2[k + 1] = otu_label[OTU2];
		edge1[k] = edge1[k + 1] = cur_nod;

		/* get the distances between all OTUs but the 2 selected ones
		   and the latter:
		   a) get the sum for both
		   b) compute the distances for the new OTU */
                double miu=0;
                double miuSum=0;
                for(i=1;i<=n;i++)
                 {
                   if(i == OTU1 || i==OTU2)continue;
                   //Rprintf("index(%i,%i)=%i index(%i,%i)=%i",i,OTU1,give_index(i,OTU1,n),i,OTU2,give_index(i,OTU2,n));
                   miuSum+=(1/(v[give_index(i,OTU1,n)]+v[give_index(i,OTU2,n)]));
                 }
                miuSum=1/miuSum;
                miu=miuSum/2;

                double eLenSum=0;
                for(i=1;i<=n;i++)
                 {
                   if(i == OTU1 || i==OTU2)continue;

                   double wi=miu/(v[give_index(i,OTU1,n)]+v[give_index(i,OTU2,n)]);
                   eLenSum+=wi*(D[give_index(i,OTU1,n)]-D[give_index(i,OTU2,n)]);
                 }

                edge_length[k]=D[give_index(OTU1,OTU2,n)]/2 + eLenSum;

                eLenSum=0;
                /*for(i=1;i<=n;i++)
                 {
                   if(i == OTU1 || i==OTU2)continue;

                   double wi=miu/(v[give_index(i,OTU1,n)]+v[give_index(i,OTU2,n)]);
                   eLenSum+=wi*(D[give_index(i,OTU2,n)]-D[give_index(i,OTU1,n)]);
		 }*/

                edge_length[k+1]=D[give_index(OTU1,OTU2,n)] - edge_length[k];

		A = D[smallest];
		ij = 0;
		for (i = 1; i <= n; i++) {
			if (i == OTU1 || i == OTU2) continue;
			double xi = D[give_index(i, OTU1, n)]; /* dist between OTU1 and i */
 			double yi = D[give_index(i, OTU2, n)]; /* dist between OTU2 and i */
                        double lamb=v[give_index(i,OTU2,n)]/(v[give_index(i,OTU2,n)]+v[give_index(i,OTU1,n)]);
			new_dist[ij] = lamb*(xi-edge_length[k])+(1-lamb)*(yi-edge_length[k+1]);
                        new_v[ij]=(v[give_index(i,OTU2,n)]*v[give_index(i,OTU1,n)])/(v[give_index(i,OTU2,n)]+v[give_index(i,OTU1,n)]);
			ij++;
		}
		/* compute the branch lengths */
                //Rprintf("l2:%f \n",edge_length[k+1]);
		/* update before the next loop
		   (we are sure that OTU1 < OTU2) */
		if (OTU1 != 1)
			for (i = OTU1; i > 1; i--)
				otu_label[i] = otu_label[i - 1];
		if (OTU2 != n)
			for (i = OTU2; i < n; i++)
				otu_label[i] = otu_label[i + 1];
		otu_label[1] = cur_nod;

		for (i = 1; i < n; i++) {
			if (i == OTU1 || i == OTU2) continue;
			for (j = i + 1; j <= n; j++) {
				if (j == OTU1 || j == OTU2) continue;
				new_dist[ij] = D[DINDEX(i, j)];
                                new_v[ij]=v[give_index(i,j,n)];
				ij++;
			}
		}

		n--;
		for (i = 0; i < n*(n - 1)/2; i++)
                 {D[i] = new_dist[i];
                  v[i] = new_v[i];
                 }
		cur_nod--;
		k = k + 2;
	}

	for (i = 0; i < 3; i++) {
		edge1[*N*2 - 4 - i] = cur_nod;
		edge2[*N*2 - 4 - i] = otu_label[i + 1];
	}

	edge_length[*N*2 - 4] = (D[0] + D[1] - D[2])/2;
	edge_length[*N*2 - 5] = (D[0] + D[2] - D[1])/2;
	edge_length[*N*2 - 6] = (D[2] + D[1] - D[0])/2;
}
