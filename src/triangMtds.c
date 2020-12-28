/* triangMtds.c    2012-04-02 */

/* Copyright 2011-2012 Andrei-Alin Popescu */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

void C_triangMtds(double* d, int* np, int* ed1,int* ed2, double* edLen)
 {
    int n=*np;
    int k=0;
    int i=0;
    int j=0;
    int ij=-1;
    int m[n+1];
    int c[n+1];
    int s[n+1];
    int o[n+1];
    for(i=1;i<=n;i++)
    {m[i]=0;
     c[i]=i;
     s[i]=0;
     for(j=1;j<=n;j++)
      {
        if(i==j){m[i]++;continue;}
        if(d[give_index(i,j,n)]==-1)continue;
        m[i]++;
      }
    }

   for(i=1;i<n;i++)
    for(j=i+1;j<=n;j++)
      {
        if(m[i]<m[j])
         {
            int aux=m[i];
            m[i]=m[j];
            m[j]=aux;
            aux=c[i];
            c[i]=c[j];
            c[j]=aux;
         }
      }
   ij=1;
   while(m[ij]==n){ij++;}
   for(i=ij;i<n;i++)
     {if(s[c[i]]==1)continue;
      for(j=i+1;j<=n;j++)
      {
       if(s[c[j]]==1)continue;
       if(d[give_index(c[i],c[j],n)]==-1)
         {s[c[j]]=1;
         }
      }
     }

   for(i=1;i<=n;i++)
     {
       if(s[i]==0)
         {
          k++;
          o[k]=i;
         }
     }
   double* sub_d=(double*)R_alloc(k*(k - 1)/2, sizeof(double));
   ij=0;
   for(i=1;i<n;i++)
   {if(s[i]==1)continue;
    for(j=i+1;j<=n;j++)
     {
      if(s[j]==1)continue;
      sub_d[ij++]=d[give_index(i,j,n)];
     }
   }

   C_triangMtd(sub_d,&k,ed1,ed2,edLen);
   for(i=0;i<2*k-3;i++)
     {
       if(ed1[i]>k)
        {
          ed1[i]+=(n-k);
        }
       if(ed2[i]>k)
        {
          ed2[i]+=(n-k);
        }
     }
   for(i=0;i<2*k-3;i++)
     {
      if(ed2[i]<=k)
       {
          ed2[i]=o[ed2[i]];
       }
     }
   for(i=1;i<=n;i++)
    {
      if(s[i]==0)continue;//take only leaves not in Y
      m[i]=0;
      for(j=1;j<=n;j++)
       {
         if(s[j]==1)continue;//take only leaves already in Y
         if(d[give_index(i,j,n)]==-1)continue;//igonore if distance unknown
         m[i]++;
       }
    }
 int numEdges=2*k-4;//0-based, so subtract 1
 //Rprintf("numEdge=%i",numEdges);
 int nv=(k-2)+n;
 while(k<n)
 {
   //now find element in X\Y such that it has most known distances to already
   //built tree until all leaves are added or we can not find a place to attach
   //the new element
   //s[i]=1 => i not added to tree
    int max=-1;
    int maxPos=-1;
    for(i=1;i<=n;i++)
    {
     if(s[i]==0)continue;
     if(m[i]>max)
      {
        max=m[i];
        maxPos=i;
      }
    }
   s[maxPos]=0;//mark maxPos as added
   //calculate new m values for leaves not added, i.e we just increment any
   //already present value by 1 if we know the distance between i and maxPos
   for(i=1;i<=n;i++)
    {
      if(s[i]==0)continue;
      if(d[give_index(i,maxPos,n)]==-1)continue;
      m[i]++;
    }

   //find path to attach maxPos to, grow tree
        double minDist=1e50;
        int z=maxPos;
        int x=-1,y=-1;
        for(i=1;i<n;i++)
        {if(s[i]==1 || d[give_index(i,z,n)]==-1 || i==z)continue;
         for(j=i+1;j<=n;j++)
          {
            if(s[j]==1 || d[give_index(j,z,n)]==-1 || j==z)continue;
            double tDist=(d[give_index(i,z,n)]+d[give_index(j,z,n)]-d[give_index(i,j,n)])/2;
            if(tDist<minDist)
             {
                minDist=tDist;
                x=i;
                y=j;
             }
          }
        }
        if(x==-1 || y==-1)
         {error("could not build tree from given distance matrix");
         }
        int* ord=getPathBetween(x,y,n,ed1,ed2,numEdges);
        /*for(i=1;i<=2*n-2;i++)
         {Rprintf("ord[%i]=%i ",i,ord[i]);
         }
        Rprintf("\n");*/
        //look for the edge on the path x to y to subdivide

        int p=x;
        double sum=0;
        double prevSum=0;
        int aux=0;

        int subdiv=-1;//index of edge to subdivide
        //error("d[y,x]=%f,d[z,x]=%f,d[z,y]=%f\n",d[give_indexx(y,x,n)],d[give_indexx(z,x,n)],d[give_indexx(z,y,n)]);
        double lx=0.5*(d[give_indexx(y,x,n)]+d[give_indexx(z,x,n)]-d[give_indexx(z,y,n)]);//distance of attachment
       // Rprintf("adding %i on the path between %i and %i at a distance from x of %f and a distance of %f from tree",z,x,y,lx,minDist);
        //point from x
        //Rprintf("Adding leaf %i, between %i and %i\n",z,x,y);
       // Rprintf("%i situated at a distance of %d from tree",z,minDist);
        int sw=0;
        //Rprintf("path between %i and %i\n",x,y);
        //int cc=0;
        while(p!=y && sum<lx)
          { //cc++;
            aux=p;
            //error("%i -> %i ",p,ord[p]);
            p=ord[p];
            prevSum=sum;
            for(i=0;i<=numEdges;i++)
              {
                if((ed1[i]==aux && ed2[i]==p)||(ed2[i]==aux && ed1[i]==p))
                  {
                    if(ed1[i]==aux && ed2[i]==p){sw=1;}
                    subdiv=i;
                    sum+=edLen[i];
                  }
              }
            //if(cc>1000)error("failed to follow path between x=%i y=%i\n",x,y);
          }


        nv++;
        //subdivide subdiv with a node labelled nv
        //length calculation

        int edd=ed2[subdiv];
        ed2[subdiv]=nv;
        edLen[subdiv]= (sw==1)?(lx-prevSum):(sum-lx);//check which 'half' of the
                                                     //path the leaf belongs to
                                                     //and updates accordingly
        //error("sum=%f, prevsum=%f\n",sum,prevSum);
        //error("lx-prevSum=%f, sum-lx=%f, minDist=%f",lx-prevSum,sum-lx,minDist);
        //Rprintf("adding %i on path %i %i, at distance %f from %i, and %f from tree\n",z,x,y,lx,x,minDist);
       // Rprintf("subdividing edge %i\n",subdiv);
        numEdges++;
        ed1[numEdges]=nv;
        ed2[numEdges]=edd;
        edLen[numEdges]= (sw==1)?(sum-lx):(lx-prevSum);
        numEdges++;
        edLen[numEdges]=minDist;
        ed1[numEdges]=nv;
        ed2[numEdges]=z;
   k++;
 }
 }
