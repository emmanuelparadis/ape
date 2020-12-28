/* additive.c    2011-10-11 */

/* Copyright 2011 Andrei-Alin Popescu */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

void C_additive(double *dd, int* np, int* mp, double *ret)//d received as dist object, -1 for missing entries
{
    int n=*np;
    int m=*mp;
    int i=0,j=0;
    double max=dd[0];
    double d[n][n];
    for(i=1;i<n;i++)
    {d[i-1][i-1]=0;
     for(j=i+1;j<=n;j++)
      {
         d[i-1][j-1]=d[j-1][i-1]=dd[give_index(i,j,n)];
         if(dd[give_index(i,j,n)]>max)
          {
            max=dd[give_index(i,j,n)];
          }
      }
    }
    d[n-1][n-1]=0;

  int entrCh=0;
   do{
    entrCh=0;
    for(i=0;i<n-1;i++)
     for(j=i+1;j<n;j++)
      {
         if(d[i][j]!=-1)continue;
         double minimax=max;
         int k=0;
         int sw=0;
         for(k=0;k<n;k++)
          {int l=0;
             if(d[i][k]==-1 || d[j][k]==-1)continue;
           for(l=0;l<n;l++)
            {
             if(k==l || d[k][l]==-1 || d[i][l]==-1 || d[j][l]==-1)continue;
             sw=1;
             double mx=(((d[i][k]+d[j][l])>(d[i][l]+d[j][k]))?(d[i][k]+d[j][l]):(d[i][l]+d[j][k]));
             mx-=d[k][l];
             if(mx<minimax){minimax=mx;}
            }
          }
        if(sw==1)
          {
            d[i][j]=d[j][i]=minimax;
            m--;
            entrCh=1;
          }
      }
   }while(entrCh==1);
  int ij=0;
  for(i=0;i<n;i++)
   for(j=0;j<n;j++)
    {
       ret[ij++]=d[i][j];
    }
}
