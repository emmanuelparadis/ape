'biplot.pcoa' <- 
    function(x, Y=NULL, plot.axes=c(1,2), dir.axis1=1, dir.axis2=1,rn=NULL,main=NULL, ...)
# x = output object from function pcoa.R
# Y = optional sites-by-variables data table
# plot.axes = the two axes to be plotted
# dir.axis.1 = -1 to revert axis 1 for the projection of points and variables
# dir.axis.2 = -1 to revert axis 2 for the projection of points and variables
# rn = an optional vector, length n, of object name labels
# Customize the title of the biplot with argument 'main'. Ex.: main="My own PCoA title".
#
# Corrected version, March 2017 - This version draws biplots from the principal coordinates (x$vectors.cor) with Lingoes or Cailliez correction, when applicable.
#
# Author: Pierre Legendre, January 2009, March 2017
	{
	if (!inherits(x, "pcoa")) stop("Object of class 'pcoa' expected")
	pr.coo <- x$vectors
	if(x$correction[2] > 1) pr.coo <- x$vectors.cor
	k <- ncol(pr.coo)
	if(k < 2) stop("There is a single eigenvalue. No plot can be produced.")
	if(k < plot.axes[1]) stop("Axis",plot.axes[1],"does not exist.")
	if(k < plot.axes[2]) stop("Axis",plot.axes[2],"does not exist.")

	if(!is.null(rn)) rownames(pr.coo) <- rn
	labels = colnames(pr.coo[,plot.axes])
	diag.dir <- diag(c(dir.axis1,dir.axis2))
	pr.coo[,plot.axes] <- pr.coo[,plot.axes] %*% diag.dir

	if(is.null(Y)) {
		limits <- apply(pr.coo[,plot.axes], 2, range) 
		ran.x <- limits[2,1] - limits[1,1]
		ran.y <- limits[2,2] - limits[1,2]
		xlim <- c((limits[1,1]-ran.x/10), (limits[2,1]+ran.x/5)) 
		ylim <- c((limits[1,2]-ran.y/10), (limits[2,2]+ran.y/10))

		par(mai = c(1.0, 1.0, 1.0, 0.5))
		plot(pr.coo[,plot.axes],xlab=labels[1],ylab=labels[2],xlim=xlim,ylim=ylim,asp=1)
		text(pr.coo[,plot.axes], labels=rownames(pr.coo), pos=4, cex=1, offset=0.5)
		if(is.null(main)) {
			title(main = "PCoA ordination", line=2)
   			} else {
   			title(main = main, family="serif", line=2)
   			}
	
		} else {
		# Find positions of variables in biplot:
		# construct U from covariance matrix between Y and standardized point vectors
		# (equivalent to PCA scaling 1, since PCoA preserves distances among objects)
		n <- nrow(Y)
		points.stand <- scale(pr.coo[,plot.axes])
		S <- cov(Y, points.stand)
		U <- S %*% diag((x$values$Eigenvalues[plot.axes]/(n-1))^(-0.5))
		colnames(U) <- colnames(pr.coo[,plot.axes])

		par(mai = c(1, 0.5, 1.4, 0))
		biplot(pr.coo[,plot.axes], U, xlab=labels[1], ylab=labels[2])
		if(is.null(main)) {
			title(main = c("PCoA biplot","Response variables projected","as in PCA with scaling 1"), line=4)
   			} else {
   			title(main = main, family="serif")
   			}
	}
    invisible()
}
