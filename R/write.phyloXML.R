## write.phyloXML.R (2024-06-17)

##   Write Tree File in PhyloXML Format

## Copyright 2002-2024 Emmanuel Paradis, Federico Marotta

## This file is part of the R-package `ape'.
## See the file ../COPYING for licensing issues.

# Save a "phylo" object to a file in phyloXML format
#
# tree: An object of class "phylo" or "multiPhylo".
write.phyloXML <- function(phy, file = "", tree.names = FALSE) {
  phyloxml <- phylo_to_xml(phy, tree.names)
  cat(as.character(phyloxml), file = file)
}

phylo_to_xml <- function(phy, tree.names = FALSE) {
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Please install the `xml2` package if you want to write phyloXML files.")
  }
  if (inherits(phy, "phylo")) {
    phy <- c(phy)
  }
  n_trees <- length(phy)
  if (is.null(attr(phy, "TipLabel"))) {
    for (i in seq_len(n_trees)) {
      phy[[i]]$tip.label <- checkLabel(phy[[i]]$tip.label)
    }
  } else {
    attr(phy, "TipLabel") <- checkLabel(attr(phy, "TipLabel"))
    phy <- .uncompressTipLabel(phy)
  }
  if (is.logical(tree.names)) {
    if (tree.names) {
      tree.names <- if (is.null(names(phy))) {
        paste0("tree", seq_len(n_trees))
      } else {
        names(phy)
      }
    } else {
      tree.names <- character(n_trees)
    }
  }
  phyloxml <- xml2::xml_new_root("phyloxml",
    `xmlns:xsi` = "http://www.w3.org/2001/XMLSchema-instance",
    xmlns = "http://www.phyloxml.org",
    `xsi:schemaLocation` = "http://www.phyloxml.org http://www.phyloxml.org/1.20/phyloxml.xsd"
  )
  lapply(seq_len(n_trees), function(i) {
    root_idx <- unique(
      phy[[i]]$edge[! phy[[i]]$edge[, 1] %in% phy[[i]]$edge[, 2], 1]
    )
    stopifnot(length(root_idx) == 1)
    clades <- .phylo_to_xml_clades(root_idx, phy[[i]])
    if (!is.null(phy[[i]]$root.edge)) {
      xml2::xml_set_attr(clades, "branch_length", phy[[i]]$root.edge)
    }
    phylogeny <- xml2::read_xml("<phylogeny></phylogeny>")
    xml2::xml_set_attr(phylogeny, "rooted", tolower(is.rooted(phy[[i]])))
    if (nchar(tree.names[i])) {
      xml2::xml_add_child(phylogeny, "name", tree.names[i])
    }
    xml2::xml_add_child(phylogeny, clades)
    xml2::xml_add_child(phyloxml, phylogeny)
  })
  return(phyloxml)
}

.phylo_to_xml_clades <- function(parent_idx, tree) {
  parent <- xml2::read_xml("<clade></clade>")
  node_name <- if (parent_idx <= length(tree$tip.label)) {
    tree$tip.label[parent_idx]
  } else if (!is.null(tree$node.label)) {
    tree$node.label[parent_idx - length(tree$tip.label)]
  } else {
    parent_idx
  }
  xml2::xml_add_child(parent, "name", node_name)
  which_children <- which(tree$edge[, 1] == parent_idx)
  lapply(which_children, function(which_child) {
    child_idx <- tree$edge[which_child, 2]
    child <- .phylo_to_xml_clades(child_idx, tree)
    if (!is.null(tree$edge.length)) {
      branch_length <- tree$edge.length[which_child]
      xml2::xml_set_attr(child, "branch_length", branch_length)
    }
    xml2::xml_add_child(parent, child)
  })
  return(parent)
}
