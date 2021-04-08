/* treePop.c    2013-09-19 */

/* Copyright 2011-2013 Andrei-Alin Popescu */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include <math.h>
#include <R.h>
#include <stdint.h>

int lsb(uint8_t * a)
{
	int i = 0;
	while (a[i] == 0) i++; /* count number of elements = 0 */

	int b = 7;
	while ((a[i] & (1 << b)) == 0) b--;

	return i*8 + (8 - b);
}

short count_bits(uint8_t n)
{
	short c; /* c accumulates the total bits set in v */
	for (c = 0; n; c++)
		n &= n - 1; /* clear the least significant bit set */
	return c;
}

uint8_t* setdiff(uint8_t* x, uint8_t *y, int nrow) //x-y
{
	int i = 0;
	uint8_t* ret = (uint8_t*)R_alloc(nrow, sizeof(uint8_t));
	for (i = 0; i < nrow; i++) {
		uint8_t tmp = (~y[i]);

		ret[i] = (x[i] & tmp);
	}
	return ret;
}

void C_treePop(int* splits, double* w,int* ncolp,int* np, int* ed1, int* ed2, double* edLen)
  {
    int n=*np;
    int ncol=*ncolp;
    int nrow=ceil(n/(double)8);
    uint8_t split[nrow][ncol];
    int i=0,j=0;

    for(i=0;i<ncol;i++)
    {
     for(j=0;j<nrow;j++)
       {
         split[j][i]=(uint8_t)splits[i*nrow+j];
       }
    }

   uint8_t vlabs[2*n-1][nrow];
   for(i=0;i<nrow-1;i++)
    {
       vlabs[n+1][i]=255;
    }
   vlabs[n+1][nrow-1]=~((uint8_t)(pow(2,8-(n%8))-1));

   int bitt_count[ncol];
   uint8_t msk=~((uint8_t)(pow(2,8-(n%8))-1));//mask out trailing bits
   for(i=0;i<ncol;i++)
    {
      int sum=0;
      for(j=0;j<nrow-1;j++)
       {
         sum+=count_bits(split[j][i]);
       }
      uint8_t bt=split[nrow-1][i];
      bt&=msk;
      sum+=count_bits(bt);
      if(sum>n/2)
       {
         for(j=0;j<nrow;j++)
          {
             split[j][i]=~split[j][i];
          }
         split[nrow-1][i]&=msk;
         sum=n-sum;
       }
      bitt_count[i]=sum;
    }
   int ind[ncol];
   for(i=0;i<ncol;i++)
     {
       ind[i]=i;
     }

   for(i=0;i<ncol-1;i++)
     {
       for(j=i+1;j<ncol;j++)
        {
          if(bitt_count[i]<bitt_count[j])
            { int aux;
              aux=bitt_count[i];
              bitt_count[i]=bitt_count[j];
              bitt_count[j]=aux;
              aux=ind[i];
              ind[i]=ind[j];
              ind[j]=aux;
            }
        }
     }
   int nNodes=n+1;
   int numEdges=0;
   for(i=0;i<ncol;i++)
    {  int ii=0;
       uint8_t sp[nrow];
       for(ii=0;ii<nrow;ii++)
         {//copy split into sp
          sp[ii]=split[ii][ind[i]];
         }
      //search for node whose labellings are a superset of the current split
      for(j=n+1;j<=nNodes;j++)
       {
          uint8_t vl[nrow];
          for(ii=0;ii<nrow;ii++)
            {//copy vertex labeling into vl
              vl[ii]=vlabs[j][ii];
            }
	  int sw=0;
         uint8_t* sd=setdiff(sp,vl,nrow);
         for(ii=0;ii<nrow/*-1*/;ii++)
           {
             if(sd[ii]!=0)sw=1;
           }

         sd[nrow-1]&=msk;
         if(sd[nrow-1]!=0)sw=1;
         if(sw==0)//setdiff==0, so we split vertex j
          {
             ed1[numEdges]=j;
             int gn=-1;
             if(bitt_count[i]>=2)//if not trivial split
             {nNodes++;
              gn=nNodes;
             }
             else
             {
               gn=lsb(sp);
             }
             ed2[numEdges]=gn;
             edLen[numEdges]=w[ind[i]];
             numEdges++;
             uint8_t* sdd=setdiff(vl,sp,nrow);
             for(ii=0;ii<nrow;ii++)//label current vertex byset difference
                                   //and new vertex by actual split
               {
                 vlabs[j][ii]=sdd[ii];
                 vlabs[gn][ii]=sp[ii];
               }
             break;
          }

       }
    }
  }
