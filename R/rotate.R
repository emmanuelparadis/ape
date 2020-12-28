## rotate.R (2014-06-05)

##   Ancestral Character Estimation

## Copyright 2007 Christoph Heibl

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

rotate <- function(phy, node, polytom = c(1,2)){
	# load DESCENDANTS function
	DESCENDANTS <- function(tree, node){
		tips <- length(tree$tip.label)
		x <- tree$edge[,2][tree$edge[,1] == node]
		while(max(x) > tips){
			x <- x[x > tips]
			for(h in 1:length(x)) tree$edge <- tree$edge[!tree$edge[,2] == x[h],]
			for(i in 1:length(x)) tree$edge[,1][tree$edge[,1] == x[i]] <- node
			x <- tree$edge[,2][tree$edge[,1] == node]
			}
		x
		}

	if (!inherits(phy, "phylo")) # is phy of class phylo?
        stop("object \"phy\" is not of class \"phylo\"")
        phy <- reorder(phy) # added by EP 2014-06-05
    nb.tips <- length(phy$tip.label) # number of tiplabels
	max.int.node <- phy$Nnode+nb.tips # number of last internal node
	nb.edges <- dim(phy$edge)[1] # number of branches
	if (length(node) == 2){ # get MRCA if tips are given for node
    	if (mode(node) == "character"){
    		if (any(!node %in% phy$tip.label)) # do tiplabels correspond
        		stop("object \"node\" contains tiplabels not present in object \"phy\"")
    		tips <- cbind(phy$tip.label, 1:nb.tips)
    		node[1] <- tips[,2][tips[,1] == node[1]]
			node[2] <- tips[,2][tips[,1] == node[2]]
			node <- as.numeric(node)
    		}
    	if (any(!node %in% 1:nb.tips))
        	stop("object \"node\" does not contain terminal nodes")
    	node <- getMRCA(phy, node)
    	}
	if (node  <= nb.tips || node > max.int.node) # is node really internal?
        stop("object \"node\" is not an internal node of object \"phy\"")
	with.br.length <- !is.null(phy$edge.length) # does phy contain brlength?
	G <- cbind(phy$edge, 1:(length(phy$edge)/2))
	N <- phy$edge[phy$edge[,1] == node]
	N <- N[N != node]
	if (length(N) > 2) N <- N[polytom]
	CLADE1 <- N[1]
	CLADE2 <- N[2]
# do clades comprise interior nodes?
	if (CLADE1 > nb.tips) CLADE11 <- DESCENDANTS(phy, CLADE1)
	if (CLADE2 > nb.tips) CLADE22 <- DESCENDANTS(phy, CLADE2)
# calculate inidices of clades in phy.edge
		if (CLADE1 > nb.tips){
			c1 <- G[,3][G[,2] == CLADE1]
			c2 <- G[,3][G[,2] == max(CLADE11)]
			} else {
			c1 <- G[,3][G[,2] == CLADE1]
			c2 <- G[,3][G[,2] == CLADE1]
			}
		if (CLADE2 > nb.tips){
			c3 <- G[,3][G[,2] == CLADE2]
			c4 <- G[,3][G[,2] == max(CLADE22)]
			} else {
			c3 <- G[,3][G[,2] == CLADE2]
			c4 <- G[,3][G[,2] == CLADE2]
			}

# create new phy$edge and  phy$edge.length
if (c2+1 == c3){
	if (c1 == 1 && c4 != nb.edges){
		phy$edge <- rbind(phy$edge[c3:c4,], phy$edge[c1:c2,], phy$edge[(c4+1):nb.edges,])
			if (with.br.length)
			phy$edge.length <- c(phy$edge.length[c3:c4], phy$edge.length[c1:c2], phy$edge.length[(c4+1):nb.edges])
		}
	if (c1 !=1 && c4 == nb.edges){
		phy$edge <- rbind(phy$edge[1:(c1-1),], phy$edge[c3:c4,], phy$edge[c1:c2,])
			if (with.br.length)
			phy$edge.length <- c(phy$edge.length[1:(c1-1)], phy$edge.length[c3:c4], phy$edge.length[c1:c2])
		}
	if (c1 !=1 && c4 != nb.edges){
		phy$edge <- rbind(phy$edge[1:(c1-1),], phy$edge[c3:c4,], phy$edge[c1:c2,], phy$edge[(c4+1):nb.edges,])
			if (with.br.length)
			phy$edge.length <- c(phy$edge.length[1:(c1-1)], phy$edge.length[c3:c4], phy$edge.length[c1:c2], phy$edge.length[(c4+1):nb.edges])
		}
	if (c1 ==1 && c4 == nb.edges){
		phy$edge <- rbind(phy$edge[c3:c4,], phy$edge[c1:c2,])
			if (with.br.length)
			phy$edge.length <- c(phy$edge.length[c3:c4], phy$edge.length[c1:c2])
		}
	}
else {
	if (c1 == 1 && c4 != nb.edges){
		phy$edge <- rbind(phy$edge[c3:c4,], phy$edge[(c2+1):(c3-1),], phy$edge[c1:c2,], phy$edge[(c4+1):nb.edges,])
			if (with.br.length)
			phy$edge.length <- c(phy$edge.length[c3:c4], phy$edge.length[(c2+1):(c3-1)], phy$edge.length[c1:c2], phy$edge.length[(c4+1):nb.edges])
		}
	if (c1 !=1 && c4 == nb.edges){
		phy$edge <- rbind(phy$edge[1:(c1-1),], phy$edge[c3:c4,], phy$edge[(c2+1):(c3-1),], phy$edge[c1:c2,])
			if (with.br.length)
			phy$edge.length <- c(phy$edge.length[1:(c1-1)], phy$edge.length[c3:c4], phy$edge.length[(c2+1):(c3-1)], phy$edge.length[c1:c2])
		}
	if (c1 !=1 && c4 != nb.edges){
		phy$edge <- rbind(phy$edge[1:(c1-1),], phy$edge[c3:c4,], phy$edge[(c2+1):(c3-1),], phy$edge[c1:c2,], phy$edge[(c4+1):nb.edges,])
			if (with.br.length)
			phy$edge.length <- c(phy$edge.length[1:(c1-1)], phy$edge.length[c3:c4], phy$edge.length[(c2+1):(c3-1)], phy$edge.length[c1:c2], phy$edge.length[(c4+1):nb.edges])
			}
	if (c1 ==1 && c4 == nb.edges){
		phy$edge <- rbind(phy$edge[c3:c4,], phy$edge[(c2+1):(c3-1),], phy$edge[c1:c2,])
			if (with.br.length)
			phy$edge.length <- c(phy$edge.length[c3:c4], phy$edge.length[(c2+1):(c3-1)], phy$edge.length[c1:c2])
		}
	}
	phy
	}
