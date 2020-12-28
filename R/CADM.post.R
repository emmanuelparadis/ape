`CADM.post` <-
	function(Dmat, nmat, n, nperm=99, make.sym=TRUE, weights=NULL, mult="holm", mantel=FALSE, silent=FALSE)
{
### Function to carry out a posteriori tests of the contribution of individual 
### matrices to the congruence of a group of distance matrices.
###
### copyleft - Pierre Legendre, December 2008
###
### Reference -
### Legendre, P. and F.-J. Lapointe. 2004. Assessing congruence among distance 
### matrices: single malt Scotch whiskies revisited. Australian and New Zealand 
### Journal of Statistics 46: 615-629.
###
### Parameters of the function --
###
### Dmat = A text file listing the distance matrices one after the other, with
###        or without blank lines.
###        Each matrix is in the form of a square distance matrix with 0's 
###        on the diagonal.
###
### nmat = number of distance matrices in file Dmat.
###
### n = number of objects in each distance matrix. All matrices have same n.
###
### nperm = number of permutations for the tests.
###
### make.sym = TRUE: turn asymmetric matrices into symmetric matrices by 
###            averaging the two triangular portions.
###          = FALSE: analyse asymmetric matrices as they are.
###
### weights = a vector of positive weights for the distance matrices. 
###           Example: weights = c(1,2,3)
###         = NULL (default): all matrices have same weight in calculation of W.
###
### mult = method for correcting P-values due to multiple testing. The methods 
###        are "holm" (default), "sidak", and "bonferroni". The Bonferroni 
###        correction is overly conservative; it is not recommended. It is 
###        included to allow comparisons with the other methods.
###
### mantel = TRUE: Mantel statistics are computed from ranked distances,
###          as well as permutational P-values.
###        = FALSE (default): Mantel statistics and tests are not computed.
###
### silent = TRUE: informative messages will not be printed, except stopping 
###          messages. Option useful for simulation work.
###        = FALSE: informative messages will be printed.
###
################################################################################
	
	mult <- match.arg(mult, c("sidak", "holm", "bonferroni"))
	if(nmat < 2) 
		stop("Analysis requested for a single D matrix: CADM is useless")
	
	a <- system.time({

    ## Check the input file
    if(ncol(Dmat) != n) 
    	stop("Error in the value of 'n' or in the D matrices themselves")
    nmat2 <- nrow(Dmat)/n
    if(nmat2 < nmat)  # OK if 'nmat' < number of matrices in the input file
    	stop("Number of input D matrices = ",nmat2,"; this value is < nmat")

    nd <- n*(n-1)/2
    if(is.null(weights)) {
    	w <- rep(1,nmat)
    	} else {
    	if(length(weights) != nmat) 
    		stop("Incorrect number of values in vector 'weights'")
    	if(length(which(weights < 0)) > 0) 
    		stop("Negative weights are not permitted")
    	w <- weights*nmat/sum(weights)
    	if(!silent) cat("Normalized weights =",w,'\n')
    	}
    
    ## Are asymmetric D matrices present?
    asy <- rep(FALSE, nmat)
	asymm <- FALSE
    end <- 0
    for(k in 1:nmat) {
        begin <- end+1
        end <- end+n
        D.temp <- Dmat[begin:end,]
        if(sum(abs(diag(as.matrix(D.temp)))) > 0) 
        	stop("Diagonal not 0: matrix #",k," is not a distance matrix")
        vec1 <- as.vector(as.dist(D.temp))
        vec2 <- as.vector(as.dist(t(D.temp)))
        if(sum(abs((vec1-vec2))) > 0) {
        	if(!silent) cat("Matrix #",k," is asymmetric",'\n')
        	asy[k] <- TRUE
        	asymm <- TRUE
        	}
        }
    D1 <- as.list(1:nmat)
    if(asymm) {
    	if(make.sym) {
    		if(!silent) cat("\nAsymmetric matrices were transformed to be symmetric",'\n')
    		} else {
    		nd <- nd*2
    		if(!silent) cat("\nAnalysis carried out on asymmetric matrices",'\n')
    		D2 <- as.list(1:nmat)
    		}
    	} else {
    	if(!silent) cat("Analysis of symmetric matrices",'\n')
    	}
    Y <- rep(NA,nd)
    
    ## String out the distance matrices (vec) and assemble them as columns into matrix 'Y'
    ## Construct also matrices of ranked distances D1[[k]] and D2[[k]] for permutation test
    end <- 0
    for(k in 1:nmat) {
        begin <- end+1
        end <- end+n
        D.temp <- as.matrix(Dmat[begin:end,])
        vec <- as.vector(as.dist(D.temp))
        if(asymm) {
        	if(!make.sym) {
        		## Analysis carried out on asymmetric matrices: 
        		## The ranks are computed on the whole matrix except the diagonal values.
        		## The two halves are stored as symmetric matrices in D1[[k]] and D2[[k]]
        		vec <- c(vec, as.vector(as.dist(t(D.temp))))
        		diag(D.temp) <- NA
                        D.temp2 <- rank(D.temp)
                        dim(D.temp2) <- dim(D.temp)   # Correction E. Paradis, 08may17
        		diag(D.temp2) <- 0
        		# cat("nrow =",nrow(D.temp2)," ncol =",ncol(D.temp2),'\n')
				# cat("Matrix ",k," min =",min(D.temp2)," max =",max(D.temp2),'\n')
				# cat("Matrix ",k," max values #",which(D.temp2 == max(D.temp2)),'\n')
        		D1[[k]] <- as.matrix(as.dist(D.temp2))
        		D2[[k]] <- as.matrix(as.dist(t(D.temp2)))
        		} else {
        		## Asymmetric matrices transformed to be symmetric, stored in D1[[k]]
        		vec <- (vec + as.vector(as.dist(t(D.temp)))) / 2
				D.temp2 <- (D.temp + t(D.temp)) / 2
				D.temp2 <- as.dist(D.temp2)
        		D.temp2[] <- rank(D.temp2)
				D.temp2 <- as.matrix(D.temp2)
        		D1[[k]] <- D.temp2
        		}
        	} else {
        	## Symmetric matrices are stored in D1[[k]]
        	D.temp2 <- as.dist(D.temp)
        	D.temp2[] <- rank(D.temp2)
        	D1[[k]] <- as.matrix(D.temp2)
        	}
        Y <- cbind(Y, vec)
        }
    Y <- as.matrix(Y[,-1])
    colnames(Y) <- colnames(Y,do.NULL = FALSE, prefix = "Dmat.")
    
    ## Begin calculations: compute reference value of S

		## Transform the distances to ranks, by column
		Rmat <- apply(Y,2,rank)

		## Compute the S = Sum-of-Squares of the row-marginal sums of ranks (eq. 1a)
		## The ranks are weighted during the sum by the vector of matrix weights 'w'
		sumRanks <- as.vector(Rmat%*%w)
		S <- (nd-1)*var(sumRanks)

    ## Begin a posteriori tests of individual matrices

    ## Statistics displayed for each matrix: "Mantel.mean" and "W.per.matrix"
	## Calculate the mean of the Mantel correlations on ranks for each matrix
	Mantel.cor <- cor(Rmat)
	diag(Mantel.cor) <- 0
	spear.mean <- as.vector(Mantel.cor%*%w)/(nmat-1)
	## Calculate Kendall's W for each variable
	## W.var <- ((nmat-1)*spear.mean+1)/nmat

	## P-value for each matrix: test of S, permuting values in matrix[[k]] only
	## as in program CADM.f (2004)
	## Initialize the counters
	counter <- rep(1,nmat)

	## Test each matrix 'k' in turn
	for(k in 1:nmat) {
		## Create a new Rmat table where the permuted column has been removed
		Rmat.mod <- Rmat[,-k]
				
		## Permutation loop: string out permuted matrix 'k' only
		for(j in 1:nperm) {
			order <- sample(n)
			if(asymm & !make.sym) {
				## For asymmetric matrices: permute the values within each triangular 
				## portion, stored as square matrices in D1[[]] and D2[[]]
				vec <- as.vector(as.dist(D1[[k]][order,order]))
			    vec <- c(vec, as.vector(as.dist(D2[[k]][order,order])))
				} else {
				vec <- as.vector(as.dist(D1[[k]][order,order]))
				}
			Rmat.perm <- cbind(Rmat.mod, vec)
			S.perm <- (nd-1)*var(as.vector(Rmat.perm%*%w))
			if(S.perm >= S) counter[k] <- counter[k]+1
		}
	}

	## Calculate P-values
	counter <- counter/(nperm+1)
	
	## Correction to P-values for multiple testing
		if(mult == "sidak") {
			vec.corr = NA
			for(i in 1:nmat) vec.corr = c(vec.corr, (1-(1-counter[i])^nmat))
			vec.corr <- vec.corr[-1]
			}
		if(mult == "holm") vec.corr <- p.adjust(counter, method="holm")
		if(mult == "bonferroni") vec.corr <- p.adjust(counter, method="bonferroni")

	## Create a data frame containing the results
		# table <- rbind(spear.mean, W.var, counter, vec.corr)
		# rownames(table) <- c("Mantel.mean", "W.per.matrix", "Prob", "Corrected prob")
		table <- rbind(spear.mean, counter, vec.corr)
		rownames(table) <- c("Mantel.mean", "Prob", "Corrected.prob")
		colnames(table) <- colnames(table,do.NULL = FALSE, prefix = "Dmat.")
	
	## Mantel tests
	if(mantel) {
		diag(Mantel.cor) <- 1
		rownames(Mantel.cor) <- colnames(table)
		colnames(Mantel.cor) <- colnames(table)
		Mantel.prob <- matrix(1,nmat,nmat)
		rownames(Mantel.prob) <- colnames(table)
		colnames(Mantel.prob) <- colnames(table)
		
		for(j in 1:nperm) {   # Each matrix is permuted independently
	    	                  # There is no need to permute the last matrix
			Rmat.perm <- rep(NA,nd)
			##
			if(asymm & !make.sym) {
				## For asymmetric matrices: permute the values within each triangular 
				## portion, stored as square matrices in D1[[]] and D2[[]]
				for(k in 1:(nmat-1)) {
					order <- sample(n)
					vec <- as.vector(as.dist(D1[[k]][order,order]))
			    	vec <- c(vec, as.vector(as.dist(D2[[k]][order,order])))
					Rmat.perm <- cbind(Rmat.perm, vec)
					}
					vec <- as.vector(as.dist(D1[[nmat]]))
			    	vec <- c(vec, as.vector(as.dist(D2[[nmat]])))
					Rmat.perm <- cbind(Rmat.perm, vec)
				} else {
				for(k in 1:(nmat-1)) {
					order <- sample(n)
					vec <- as.vector(as.dist(D1[[k]][order,order]))
					Rmat.perm <- cbind(Rmat.perm, vec)
					}
					vec <- as.vector(as.dist(D1[[nmat]]))
					Rmat.perm <- cbind(Rmat.perm, vec)
				}
			# Remove the first column of Rmat.perm containing NA
			Rmat.perm <- as.matrix(Rmat.perm[,-1])
			# Compute Mantel correlations on ranks under permutation
			Mantel.cor.perm <- cor(Rmat.perm)
			for(j2 in 1:(nmat-1)) { # Compute prob in the upper tail
				for(j1 in (j2+1):nmat) {
					if(Mantel.cor.perm[j1,j2] >= Mantel.cor[j1,j2]) Mantel.prob[j1,j2] <- Mantel.prob[j1,j2]+1
					}
				}
			}
		Mantel.prob <- as.matrix(as.dist(Mantel.prob/(nperm+1)))
		diag(Mantel.prob) <- NA # Corrected 08feb13
		}
	
	})
	a[3] <- sprintf("%2f",a[3])
	if(!silent) cat("Time to compute a posteriori tests (per matrix) =",a[3]," sec",'\n')

	out <- list(A_posteriori_tests=table, Correction.type=mult)

	if(mantel) {
		out$Mantel.cor  <- Mantel.cor
		out$Mantel.prob <- Mantel.prob
		}
	out$nperm <- nperm
	class(out) <- "CADM.post"
	out
}
