library(ape)
library(tidyverse)
library(assertthat)
library(phangorn)

# This script defines two functions and at the end demonstrate them with some example data. The goal is to reduce a tree to a set of tips specified in a list. However, this list includes internal nodes and tips, so ape::drop/keep.tip can't be used.

# The first function, "fun_Descendants_node_labels", is a wrapped for phangorn::Descendants(). The change here is that the original function only takes node numbers for the argument node. It is often desirable to instead be able to perform this kind of operation but with node labels. This new function takes in node labels, looks up their node numbers, passes that to phangorn::Descendants(), takes the output (which is tip numbers) and looks up the labels for them. This functions outputs a flat vector, not a nested set of lists like phangorn::Descendants().

# The second function is a wrapper of kind for ape::drop/keep.tip. The original ape functions takes a list of tips and either drops them or drops all but them (depending on whether it is drop.tip or keep.tip). These functions are not able to handle node labels or positions, only tips. Sometimes it is desirable to be able to also specify internal nodes and prune away all of its descendant and make it a tip. This new function is able to handle such scenarios. You can feed this function a vector of tip and node labels, and when it encounters a node it goes through and iteratively prunes away all of its descendants that are tips until that node becomes a tip. It also checks if the vector it got makes sense in terms of none of the labels being daughters of other labels in the list.

# The second function needs the first function to work.

## DISCLAIMER
# These functions were written to solve a specific problem. The functions are written perhaps a bit unusually, with tidyverse pipes and data.frames where more "elegant" base-solutions are available. This function was written for specific use cases I have encountered. If you find them useful and would like them to be improved in terms of functionality, elegance and/or speed let me know. For another take on this problem, see Erich Round's solution here: https://github.com/erichround/glottoTrees/blob/99338ff3e62fa3f282eb19ed55a60c564fd843bc/R/topology.R#L326

#######FUNCTIONS

#### FUN 1
# wrapper function for phangorn::Descendants() such that it can take node labels instead of node positions.
fun_Descendants_node_labels <- function(tree, node, type = "tips"){
  
  #making a table with the node number of each node.label and also it's associated children
  node_labels_in_edge <- tree$node.label[tree$edge[,1]-Ntip(tree)]
  tips_nodes <- tree$edge[,2]
  
  select.tip.or.node <- function(element, tree) {
    ifelse(element < Ntip(tree)+1, tree$tip.label[element], tree$node.label[element-Ntip(tree)])
  }
  
  edge_table <- data.frame(
    "parent" = tree$edge[,1], #node number of the parent
    "parent.name" = sapply(tree$edge[,1], select.tip.or.node, tree = tree),
    "child" = tree$edge[,2], #node number of the child
    "child.name" = sapply(tree$edge[,2], select.tip.or.node, tree = tree)
  )
  
  nodes_to_keep_number <- node %>% 
    as.data.frame() %>% 
    rename(parent.name = ".") %>% 
    inner_join(edge_table, by = "parent.name") %>% 
    dplyr::select("parent") %>% 
    .[,1] %>% 
    unique()
  
  descendants_node_numbers <- phangorn::Descendants(tree, node = nodes_to_keep_number, type = type)  %>%  
    unlist() %>% 
    as.numeric()
  
  descendants_tips_labels <- descendants_node_numbers %>% 
    as.data.frame() %>% 
    rename(child = ".") %>% 
    inner_join(edge_table, by = "child") %>% 
    dplyr::select("child.name") %>% 
    .[,1] %>% 
    unique()
  
  descendants_tips_labels
}

#### FUN 2

# modification to keep.tip such that you can also specify internal nodes. All descendants of internal nodes are dropped and the internal node becomes a tip.
keep.as.tip<- function(tree, tips_and_nodes_to_keep){
  
  #for bug running example
  #  tree <- clade_tree
  
  #making a vector of the items in tips_and_nodes_to_keep that are indeed nodes
  nodes_df <- tree$node.label %>% 
    as.data.frame()   
  
  nodes_to_keep <-  nodes_df %>% 
    rename(tips_and_nodes_to_keep = ".") %>% 
    inner_join(as.data.frame(tips_and_nodes_to_keep), by = "tips_and_nodes_to_keep") %>% 
    .[,1]
  
  # Sanity checking that the list supplied of nodes and tips to keep doesn't contain daughters of each other.
  all_descendants_vec <- fun_Descendants_node_labels(tree, nodes_to_keep, type = "all")  
  
  n_descdendants_in_keep_list <- all_descendants_vec %in%   tips_and_nodes_to_keep %>% sum()
  
  if (n_descdendants_in_keep_list > 0) {
    warning("The list of nodes and tips to keep supplied contains nodes and/or tips that are desdencants of other nodes in the list. This isn't allowable. Please revise list.")
    stop()
  } else{
    
    # going through the internal nodes one by one and pruning away descendant tips. It has to be tips that are pruned each time because drop.tip() only takes tips. So the while-loop works "upwards" in the tree until it hits the node in question and then stops.
    for (node_name in nodes_to_keep) {
      #node_name <- nodes_to_keep[2]
      tips_to_drop <- fun_Descendants_node_labels(tree, node = node_name, type = "tips") 
      
      while (length(tips_to_drop) !=0){
        
        tree <- ape::drop.tip(tree, tip = tips_to_drop, trim.internal = F, collapse.singles = F) #it is very important that collapse.singles is set to F, otherwise you may lose the node you intend to keep.
        tips_to_drop <- fun_Descendants_node_labels(tree, node = node_name, type = "tips") 
      }
    }
    
    tips_to_drop <- setdiff(tree$tip.label , tips_and_nodes_to_keep)
    
    tree <- drop.tip(tree, tip = tips_to_drop, collapse.singles = F)
    
    #double checking that we are indeed keeping what we intended to keep, not more or less
    n_to_keep <- tips_and_nodes_to_keep %>% length()
    
    joined <- tree$tip.label %>% 
      as.data.frame() %>%
      rename(tips_and_nodes_to_keep = ".") %>% 
      inner_join(as.data.frame(tips_and_nodes_to_keep), by = "tips_and_nodes_to_keep")
    
    x <- assertthat::assert_that(nrow(joined) == n_to_keep, msg = "The resulting tree does not match the list of the nodes and tips to keep.")
    
    #outputting the trimmed tree
    tree
  }
}