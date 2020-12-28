`CADM.global` <-
	function(Dmat, nmat, n, nperm=99, make.sym=TRUE, weights=NULL, silent=FALSE)
{
### Function to test the overall significance of the congruence among 
### a group of distance matrices using Kendall's coefficient of concordance W.
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
### silent = TRUE: informative messages will not be printed, except stopping 
###          messages. Option useful for simulation work.
###        = FALSE: informative messages will be printed.
###
################################################################################
	
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
    
    ## Begin calculations for global test
    
    ## Compute the reference values of the statistics: W and Chi2
	## Transform the distances to ranks, by column
	Rmat <- apply(Y,2,rank)

	## Correction factors for tied ranks (eq. 3.3)
	t.ranks <- apply(Rmat, 2, function(x) summary(as.factor(x), maxsum=nd))
	TT <- sum(unlist(lapply(t.ranks, function(x) sum((x^3)-x))))
	# if(!silent) cat("TT = ",TT,'\n')
	
	## Compute the S = Sum-of-Squares of the row-marginal sums of ranks (eq. 1a)
	## The ranks are weighted during the sum by the vector of matrix weights 'w'
	## Eq. 1b cannot be used with weights; see formula for W below
	sumRanks <- as.vector(Rmat%*%w)
	S <- (nd-1)*var(sumRanks)
	
	## Compute Kendall's W (eq. 2a)
	## Eq. 2b cannot be used with weights 
	## because the sum of all ranks is not equal to m*n*(n+1)/2 in that case
	W <- (12*S)/(((nmat^2)*((nd^3)-nd))-(nmat*TT))
		
	## Calculate Friedman's Chi-square (Kendall W paper, 2005, eq. 3.4)
	Chi2 <- nmat*(nd-1)*W

	## Test the Chi2 statistic by permutation
	counter <- 1
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
		# The test is based on the comparison of S and S.perm instead of the comparison of 
		# Chi2 and Chi2.perm: it is faster that way. 
		# S, W, and Chi2 are equivalent statistics for permutation tests.
		Rmat.perm <- as.matrix(Rmat.perm[,-1])
		S.perm <- (nd-1)*var(as.vector(Rmat.perm%*%w))
		if(S.perm >= S) counter <- counter+1
	}
	prob.perm.gr <- counter/(nperm+1)

	table <- rbind(W, Chi2, prob.perm.gr)
	colnames(table) <- "Statistics"
	rownames(table) <- c("W", "Chi2", "Prob.perm")
	})
	a[3] <- sprintf("%2f",a[3])
	if(!silent) cat("\nTime to compute global test =",a[3]," sec",'\n')
#
	# if(asymm & !make.sym) { out <- list(congruence_analysis=table, D1=D1, D2=D2)
	# } else {
	   out <- list(congruence_analysis=table)
	#    }
#
	out$nperm <- nperm
	class(out) <- "CADM.global"
	out
}
