## all.equal.phylo.R (2009-07-05)
##
##     Global Comparison of two Phylogenies

## Copyright 2006 Benoit Durand
##    modified by EP for the new coding of "phylo" (2006-10-04)

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

## Recherche de la correspondance entre deux arbres
## Parcours en profondeur et en parallele des deux arbres (current et target)
## current, target: les deux arbres a comparer
## use.edge.length: faut-il comparer les longueurs de branches ?
## use.tip.label: faut-il comparer les etiquettes de feuilles ou seulement la
##	topologie des deux arbres ?
## index.return: si TRUE, retourner la matrice de correspondance entre noeuds
##	et feuilles, une matrice a deux colonnes (current et target) avec pour
##	chaque ligne des paires d'identifiants de noeuds/feuilles, tels qu'ils
##	apparaissent dans l'attribut 'edge' des objets phylo
## tolerance, scale: parametres de comparaison des longueurs de branches
##	(voir 'all.equal')
all.equal.phylo <- function(target, current,
                        use.edge.length = TRUE,
                        use.tip.label = TRUE,
                        index.return = FALSE,
                        tolerance = .Machine$double.eps ^ 0.5,
                        scale = NULL, ...) {

	same.node <- function(i, j) {
		# Comparaison de un noeud et une feuille
		if (xor(i > Ntip1, j > Ntip2)) return(NULL)
		# Comparaison de deux feuilles
		if (i <= Ntip1) {
			if (!use.tip.label) return(c(i, j))
			if (current$tip.label[i] == target$tip.label[j])
				return(c(i, j))
			return(NULL)
  		}
  		# Comparaison de deux noeuds
		i.children <- which(current$edge[, 1] == i)
		j.children <- which(target$edge[, 1] == j)
		if (length(i.children) != length(j.children)) return(NULL)
		correspondance <- NULL
		for (i.child in i.children) {
			corresp <- NULL
			for (j.child in j.children) {
				if (!use.edge.length ||
                                    isTRUE(all.equal(current$edge.length[i.child],
                                                     target$edge.length[j.child],
                                                     tolerance = tolerance,
                                                     scale = scale)))
                                    corresp <- same.node(current$edge[i.child, 2],
                                                         target$edge[j.child, 2])
				if (!is.null(corresp)) break
			}
			if (is.null(corresp)) return(NULL)
			correspondance <- c(correspondance, i, j, corresp)
			j.children <- j.children[j.children != j.child]
		}
		return(correspondance)
	}

        Ntip1 <- length(target$tip.label)
        Ntip2 <- length(current$tip.label)
        root1 <- Ntip1 + 1
        root2 <- Ntip2 + 1
        if (root1 != root2) return(FALSE)
        ## Fix by EP so that unrooted trees are correctly compared:
        if (!is.rooted(target) && !is.rooted(current)) {
            outg <- target$tip.label[1]
            if (! outg %in% current$tip.label) return(FALSE)
            target <- root(target, outg)
            current <- root(current, outg)
        }
        ## End
	result <- same.node(root1, root2)
	if (!isTRUE(index.return)) return(!is.null(result))
	if (is.null(result)) return(result)
	result <- t(matrix(result, nrow = 2))
      colnames(result) = c('current', 'target')
      return(result)
}
