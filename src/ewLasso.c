/* ewLasso.c    2013-03-30 */

/* Copyright 2013 Andrei-Alin Popescu */

/* This file is part of the R-package `ape'. */
/* See the file ../COPYING for licensing issues. */

#include "ape.h"

int isTripletCover(int nmb, int n, int** s, int stat, int sSoFar[n], int* a)//number of sides, number of leaves, sides, side under consideration, set so far, lasso
{
	int ret=0;
	if(stat==nmb)return 1;
	int i=0;

	for(i=1;i<=n;i++)
	 {
        if(!s[stat][i])continue;//not in set
		int sw=1, j;

		for(j=1;j<=n;j++)//check if all distances to previous candidates are present
		 {
			 if(!sSoFar[j])continue;//not in set so far
			 if(!a[i*(n+1)+j]){//if not, then i is not a good candidate for this side
				   //Rprintf("failed to find distance between %i and %i, a[%i][%i]=%i \n",i,j,i,j,a[i*(n+1)+j]);
				   sw=0;
			       }
		 }

		if(!sw)continue;//not all required distances are present

		sSoFar[i]=1;//try choosing i as representative for the side
		ret+=isTripletCover(nmb,n,s,stat+1,sSoFar,a)>0?1:0;//see if, with i chosen, we can find leaves in other sides to satisfy the triplet cover condition
		sSoFar[i]=0;
	 }
  return ret;
}

void C_ewLasso(double *D, int *N, int *e1, int *e2)
{
	int n, i, j, k;

	n=*N;
	int tCov=1;
	int* a = (int*)R_alloc((n+1)*(n+1), sizeof(int));//adjacency matrix of G_{\cL} graph

	for(i=1;i<=n;i++)
	 {
	  for(j=1;j<=n;j++)
	   {
        if(D[give_index(i,j,n)]==-1)//if missing value then no edge between pair of taxa (i,j) in G
		 {
		  a[i*(n+1)+j]=a[j*(n+1)+i]=0;
		 }
		else
		 {
		  a[i*(n+1)+j]=a[j*(n+1)+i]=1;// otherwise edge between pair of taxa (i,j) in G
		 }
	   }
	 }
   //check for connectedness of G

   int *q = (int*)R_alloc(2*n-1, sizeof(int));//BFS queue
   int *v = (int*)R_alloc(2*n-1, sizeof(int));//visited?

   int p=0,u=1;//p-head of queue, u- position after last loaded element

   for(i=1;i<=n;i++)v[i]=-1;

   int stNBipartite=1, con=1, comp=1;
   int ini=1;

   /*for(i=1;i<=n;i++)
    {
	 for(j=1;j<=n;j++)
	  {
		Rprintf("a[%i][%i]=%i ",i,j,a[i*(n+1)+j]);
	  }
	 Rprintf("\n");
    }*/

   while(comp)
   {
   q[p]=ini;
   v[ini]=1;
   comp=0;
   int stNBipartiteLoc=0;//check if current connected component is bipartite
   while(p<u)//BFS
    {
     int head=q[p];
	 //Rprintf("head: %i\n",head);
	 //Rprintf("unvisited neighbours: \n");
	 for(i=1;i<=n;i++)
	  {
		if(i==head)continue;
        if(a[i*(n+1)+head]==0)continue;
		if(v[i]==v[head])//same color as vertex from which we visit-> not bipartite
		  {
			  stNBipartiteLoc=1;
		  }
		if(v[i]!=-1)continue;
		//Rprintf("vertex %i \n",i);
		q[u++]=i;
		v[i]=1-v[head];
	  }
	 p++;
    }
    stNBipartite*=stNBipartiteLoc;//anding strngly-non-bipartite over all connected components



	//check if all vertices have been visited
	for(int i=1;i<=n;i++)
	 {
	  if(v[i]==-1)
	   {
		   comp=1;
		   p=0;
		   u=1;
		   ini=i;
		   con=0;
		   break;
	   }
	 }
   }

   Rprintf("connected: %i\n",con);
   Rprintf("strongly non-bipartite: %i\n",stNBipartite);

   //finally check if \cL is triplet cover of T

   //adjencency matrix of tree, 1 to n are leaves

   int *at= (int*)R_alloc((2*n-1)*(2*n-1), sizeof(int));

   for(i=1;i<=2*n-2;i++)
    {
	  for(j=1;j<=2*n-2;j++)at[i*(2*n-1)+j]=0;
    }

   for(i=0;i<2*n-3;i++)
    {
		//Rprintf("e1[%i]=%i e2[%i]=%i \n",i,e1[i],i,e2[i]);
		at[e1[i]*(2*n-1)+e2[i]]=at[e2[i]*(2*n-1)+e1[i]]=1;
    }

  /*for(i=1;i<2*n-1;i++)
    {
	 for(j=1;j<2*n-1;j++)
	  {
		Rprintf("at[%i][%i]=%i ",i,j,at[i*(2*n-1)+j]);
	  }
	 Rprintf("\n");
    }*/

   for(i=n+1;i<=2*n-2;i++)//for each interior vertex
    {
     for(j=1;j<2*n-1;j++)//reset queue and visited veectors
       {
		v[j]=-1;
		q[j]=0;
       }

	 v[i]=1;//'disconnect' graph at i
	 int *l=(int*)R_alloc(2*n-2, sizeof(int));//vertices adjacent to i

     int nmb=0;//number of found adjacent vertices of i

	 for(j=1;j<=2*n-2;j++)//find adjacent vertices
		 {
	       if(at[i*(2*n-1)+j]==1)
		    {
			  l[nmb++]=j;
		    }
		 }

	 int** s=(int**)R_alloc(nmb,sizeof(int*));//set of leaves in each side, stored as presence/absence

	 for(j=0;j<nmb;j++)s[j]=(int*)R_alloc(n+1,sizeof(int));

	 for(j=0;j<nmb;j++)
	  {
       for(k=1;k<=n;k++)s[j][k]=0;
	  }

	 /*Rprintf("for %i \n",i);

	 for(j=0;j<nmb;j++)Rprintf("l[%i]= %i ",j,l[j]);

	 Rprintf("\n");*/

	 for(j=0;j<nmb;j++)//for each adjancet vertex, find its leafset
	   {
		 p=0;
		 u=1;
		 ini=l[j];
		 q[p]=ini;
		 v[ini]=1;
		 if(ini<=n)s[j][ini]=1;
         comp=0;
		 int nbr=0;
         while(p<u)//BFS
			{
				int head=q[p];
				//Rprintf("head: %i\n",head);
				//Rprintf("unvisited neighbours: \n");
			for(nbr=1;nbr<=2*n-1;nbr++)
			  {
				if(nbr==head)continue;
				if(at[nbr*(2*n-1)+head]==0)continue;
				if(v[nbr]!=-1)continue;
				//Rprintf("vertex %i \n",nbr);
				if(nbr<=n)s[j][nbr]=1;//if leaf, then set visited to true
				q[u++]=nbr;
				v[nbr]=1;
			  }
			p++;
		  }
	   }

	 /*Rprintf("sides for %i\n",i);
	 for(j=0;j<nmb;j++)
	  {
       int ii;
	   for(ii=1;ii<=n;ii++)
	     {
			if(s[j][ii])Rprintf("%i ",ii);
	     }
	   Rprintf("\n");
	  }*/

	 int* sSoFar= (int*)R_alloc(n+1,sizeof(int));

	 for(j=1;j<=n;j++)sSoFar[j]=0;

	 tCov*=isTripletCover(nmb,n,s,0,sSoFar,a)>0?1:0;
    }
   Rprintf("is triplet cover? %i \n",tCov);
}
