## subtreeplot.R (2017-05-26)

##  Zoom on a Portion of a Phylogeny by Successive Clicks

## Copyright 2008 Damien de Vienne

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

subtreeplot<-function(x, wait=FALSE, ...) {

    sub<-subtrees(x, wait=wait)
    y<-NULL
    plot.default(0, type="n",axes=FALSE, ann=FALSE)
    repeat {
	  split.screen(c(1,2))
        screen(2)
        if (is.null(y)) plot(x,...)
        else plot(y,sub=paste("Node :", click),...)
        screen(1)
        plot(x,sub="Complete tree",main="Type ESC or right click to exit", cex.main=0.9, ...)

        N.tip<-Ntip(x)
        N.node<-Nnode(x)

        # 5/24/17 changed by Klaus
        # coor<-plotPhyloCoor(x)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
        tips<-x$tip.label
        nodes<-x$node.label
        if (is.null(x$node.label)) nodes<-(N.tip+1):(N.tip+N.node)
        labs<-c(rep("",N.tip), nodes)

        #click<-identify(coor[,1], coor[,2], labels=labs, n=1)
        click<-identify(lastPP$xx, lastPP$yy, labels=labs, n=1)
        if (length(click) == 0) {return(y)}
        if (click > N.tip) {
            close.screen(c(1,2),all.screens = TRUE)
            split.screen(c(1,2))
            screen(1) #selects the screen to plot in
            plot(x, sub="Complete tree", ...) # plots x in screen 1 (left)
            screen(2)
            for (i in 1:length(sub)) if (sub[[i]]$name==click) break
            y<-sub[[i]]
	  }
        else cat("this is a tip, you have to choose a node\n")

      }
    on.exit(return(y))
}
