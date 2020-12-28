## read.caic.R (2005-09-21)

##   Read Tree File in CAIC Format

## Copyright 2005 Julien Dutheil

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

read.caic <- function(file, brlen=NULL, skip = 0, comment.char="#", ...)
{
  text <- scan(file = file, what = character(), sep="\n", skip = skip, comment.char = comment.char, ...)

	# Parse the whole file:
  n <- length(text) / 2
  nodes <- 1:n;
	leaf.names <- character(n)
	patterns   <- character(n)
	lengths    <- numeric(n)
	for(i in 1:n)
	{
		leaf.names[i] <- text[2*i]
		patterns[i]   <- text[2*i-1]
		lengths[i]    <- nchar(patterns[i])
	}
	# Sort all patterns if not done:
	i <- order(patterns);
	leaf.names <- leaf.names[i]
	patterns   <- patterns[i]
	lengths    <- lengths[i]

	# This inner function compares two patterns:
	test.patterns <- function(p1, p2)
	{
		t1 <- strsplit(p1, split="")[[1]]
		t2 <- strsplit(p2, split="")[[1]]
		if(length(t1) == length(t2))
		{
			l <- length(t1)
			if(l==1) return(TRUE)
			return(all(t1[1:(l-1)]==t2[1:(l-1)]) & t1[l] != t2[l])
		}
		return(FALSE)
	}

	# The main loop:
	while(length(nodes) > 1)
	{
		# Recompute indexes:
		index <- logical(length(nodes))
		maxi  <- max(lengths)
		for(i in 1:length(nodes))
		{
			index[i] <- lengths[i] == maxi
		}
		i <- 1
		while(i <= length(nodes))
		{
			if(index[i])
			{
				p <- paste("(",nodes[i],sep="")
				c <- i+1
				while(c <= length(nodes) && index[c] && test.patterns(patterns[i], patterns[c]))
				{
					p <- paste(p, nodes[c], sep=",")
					c <- c+1
				}
				if(c-i < 2) stop("Unvalid format.")
				p <- paste(p, ")", sep="")
				nodes[i]   <- p
				patterns[i]<- substr(patterns[i],1,nchar(patterns[i])-1)
				lengths[i] <- lengths[i]-1
				nodes      <- nodes   [-((i+1):(c-1))]
				lengths    <- lengths [-((i+1):(c-1))]
				patterns   <- patterns[-((i+1):(c-1))]
				index      <- index   [-((i+1):(c-1))]
			}
			i <- i+1
		}
	}

	# Create a 'phylo' object and return it:
	phy <- read.tree(text=paste(nodes[1],";", sep=""))
	phy$tip.label <- leaf.names;
	if(!is.null(brlen))
	{
		br <- read.table(file=brlen)
		phy$edge.length <- br[,1]
	}
	return(phy)
}
