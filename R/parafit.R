'parafit' <-
	function(host.D, para.D, HP, nperm=999, test.links=FALSE, seed=NULL, correction="none", silent=FALSE)
#
# Test of host-parasite coevolution
# host.D = host distance or patristic matrix (class dist or matrix)
# para.D = parasite distance or patristic matrix (class dist or matrix)
# HP = host-parasite link matrix (n.host, n.para)
#
#      Pierre Legendre, May 2009
{
epsilon <- sqrt(.Machine$double.eps)
if(is.null(seed)) {
	runif(1)
	seed <- .Random.seed[trunc(runif(1,1,626))]
	}
HP <- as.matrix(HP)

host.D <- as.matrix(host.D)
host.pc <- pcoa(host.D, correction=correction)
if(host.pc$correction[2] == 1) {
	if(min(host.pc$values[,2]) < -epsilon) stop('Host D matrix has negative eigenvalues. Rerun with correction="lingoes" or correction="cailliez"')
	sum.host.values.sq <- sum(host.pc$values[,1]^2)
	host.vectors <- host.pc$vectors
	} else {
	sum.host.values.sq <- sum(host.pc$values[,2]^2)
	host.vectors <- host.pc$vectors.cor
	}
n.host <- nrow(host.D)

para.D <- as.matrix(para.D)
para.pc <- pcoa(para.D, correction=correction)
if(para.pc$correction[2] == 1) {
	if(min(para.pc$values[,2]) < -epsilon) stop('Parasite D matrix has negative eigenvalues. Rerun with correction="lingoes" or correction="cailliez"')
	sum.para.values.sq <- sum(para.pc$values[,1]^2)
	para.vectors <- para.pc$vectors
	} else {
	sum.para.values.sq <- sum(para.pc$values[,2]^2)
	para.vectors <- para.pc$vectors.cor
	}
n.para <- nrow(para.D)

if(!silent) cat("n.hosts =", n.host, ", n.parasites =", n.para,'\n')

a <- system.time({
tracemax <- max(sum.host.values.sq, sum.para.values.sq)

if(n.host == n.para) {
	if(!silent) cat("The function cannot check if matrix HP has been entered in the right way.",'\n')
	if(!silent) cat("It will assume that the rows of HP are the hosts.",'\n')
	} else {
	temp <- dim(HP)
	if(temp[1] == n.host) {
		if(temp[2] != n.para) stop("Matrices host.D, para.D and HP not comformable")
		} else if(temp[2] == n.host) {
			if(temp[1] != n.para) stop("Matrices host.D, para.D and HP not comformable")
			HP <- t(HP)
			if(!silent) cat("Matrix HP has been transposed for comformity with host.D and para.D.",'\n')
		} else {
		stop("Matrices host.D, para.D and HP not comformable")
		}
	}
p.per.h <- apply(HP, 1, sum)
h.per.p <- apply(HP, 2, sum)
#
# Compute and test the global statistics
mat.4 <- t(host.vectors) %*% HP %*% para.vectors
global <- sum(mat.4^2)
if(nperm > 0) {
	set.seed(seed)
	nGT <- 1
	global.perm <- NA
	for(i in 1:nperm) {
		HP.perm <- apply(HP, 2, sample)
		mat.4.perm <- t(host.vectors) %*% HP.perm %*% para.vectors
		global.perm <- c(global.perm, sum(mat.4.perm^2))
		if(global.perm[i+1] >= global) nGT <- nGT+1
		}
	global.perm <- global.perm[-1]
	p.global <- nGT/(nperm+1)
	} else { p.global <- NA }

#
# Test individual H-P links
if(test.links) {
	# 1. Create the list of H-P pairs
	list.hp <- which( t(cbind(HP,rep(0,n.host))) > 0)
	HP.list <- cbind((list.hp %/% (n.para+1))+1, list.hp %% (n.para+1))
	colnames(HP.list) <- c("Host","Parasite")
	n.links <- length(list.hp)

	stat1 <- NA
	stat2 <- NA
	p.stat1 <- NA
	p.stat2 <- NA
	for(k in 1:n.links) {
		#
		# 2. Compute reference values of link statistics
		HP.k <- HP
		HP.k[HP.list[k,1], HP.list[k,2]] <- 0
		mat.4.k <- t(host.vectors) %*% HP.k %*% para.vectors
		trace.k <- sum(mat.4.k^2)
		stat1 <- c(stat1, (global-trace.k))
		den <- tracemax-global
		if(den > epsilon) { 
			stat2 <- c(stat2, stat1[k+1]/den) 
			} else { 
			stat2 <- c(stat2, NA) 
			}
		#
		# 3. Test link statistics by permutations
		if(nperm > 0) {
			set.seed(seed)
			nGT1 <- 1
			nGT2 <- 1
			nperm2 <- nperm
			#
			for(i in 1:nperm) {
				HP.k.perm <- apply(HP.k, 2, sample)
				mat.4.k.perm <- t(host.vectors) %*% HP.k.perm %*% para.vectors
				trace.k.perm <- sum(mat.4.k.perm^2)
				stat1.perm <- global.perm[i]-trace.k.perm
				if(stat1.perm >= stat1[k+1]) nGT1 <- nGT1+1
				#
				if(!is.na(stat2[k+1])) {
					den <- tracemax-global.perm[i]
					if(den > epsilon) {
						stat2.perm <- stat1.perm/den 
						if(stat2.perm >= stat2[k+1]) nGT2 <- nGT2+1
						} else {
						nperm2 <- nperm2-1
						# if(!silent) cat("In permutation #",i,"den < epsilon",'\n')
						}
					}
				}
			p.stat1 <- c(p.stat1, nGT1/(nperm+1))
			if(!is.na(stat2[k+1])) {
				p.stat2 <- c(p.stat2, nGT2/(nperm2+1))
				} else {
				p.stat2 <- c(p.stat2, NA) ### Error in previous version, corrected here
				}
			} else { 
			p.stat1 <- c(p.stat1, NA)     ### Error in previous version, corrected here
			p.stat2 <- c(p.stat2, NA)     ### Error in previous version, corrected here
			}
		}
	#
	link.table <- cbind(HP.list, stat1[-1], p.stat1[-1], stat2[-1], p.stat2[-1])
	colnames(link.table) = c("Host","Parasite","F1.stat","p.F1","F2.stat","p.F2")
	out <-list(ParaFitGlobal=global, p.global=p.global, link.table=link.table, para.per.host=p.per.h, host.per.para=h.per.p, nperm=nperm)
	} else {
	if(!silent) cat("Rerun the program with option 'test.links=TRUE' to test the individual H-P links",'\n')
	out <-list(ParaFitGlobal=global, p.global=p.global, para.per.host=p.per.h, host.per.para=h.per.p, nperm=nperm)
	}
#
})
a[3] <- sprintf("%2f",a[3])
if(!silent) cat("Computation time =",a[3]," sec",'\n')
#
class(out) <- "parafit"
out
}
