## subtrees.R (2008-04-14)

##  All subtrees of a Phylogenetic Tree

## Copyright 2008 Damien de Vienne

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

subtrees<-function(tree, wait = FALSE)
{
N.tip<-Ntip (tree)
N.node<-Nnode(tree)
limit<-N.tip+N.node
sub<-list(N.node)
u<-0

  for (k in (N.tip+1):limit)
  {
 u<-u+1
	if (wait==TRUE) cat("wait... Node",u,"out of", N.node, "treated\n")
  fils<-NULL
  pere<-res <- k
	repeat
	{
	for (i in 1: length(pere)) fils<-c(fils, tree$edge[,2][tree$edge[,1]==pere[i]])
	res<-c(res, fils)
      pere<-fils
	fils<-NULL
	if (length(pere)==0) break
	}

  len<-res[res>N.tip]
   if (u==1) {
	tree2<-tree
	len<-(N.tip+1):limit
	}
   else {
  len.tip<-res[res<N.tip+1]
  vec<-1:length(tree$tip.label)
  len.tip.stay<-setdiff(vec, len.tip)
  tree2<-drop.tip(tree, len.tip.stay)
	  }
  sub[[u]]<-tree2
  sub[[u]]$name<-k
  #sub[[u]]$name<-tree2$node.label[1]
  sub[[u]]$Ntip<-Ntip(tree2)
  sub[[u]]$Nnode<-Nnode(tree2)
  if (is.null(tree$node.label))
  	sub[[u]]$node.label<-len

  }
return(sub)
cat("\n")
}
