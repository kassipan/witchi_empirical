#!/usr/bin/env Rscript

#Function to rename labels in tree if identical
make_unique_labels <- function(labels){
  data.frame(label = labels) %>%
    group_by(label) %>%
    mutate(new_label = if (n() == 1) label else paste0(label, "_par", row_number())) %>%
    pull(new_label)
}

#Function to get offspring labels of a node 
getOffspringLabels <- function(tree, node) {
  
  if (is.na(node) || node > (length(tree$tip.label) + tree$Nnode)) {
    stop("Invalid node: ", node)
  }
  # Get the offspring of the specified node using ape's offspring function
  offspring_labels <- offspring(tree, node)
  
  # If no offspring, return an empty tibble
  if (length(offspring_labels) == 0) {
    return(tibble(label = character(0)))  # Return empty tibble if no offspring
  }
  
  # Return the labels of the offspring nodes
  offspring_labels <- tibble(label = tree$tip.label[offspring_labels]) %>%
    drop_na()
  return(offspring_labels)
}

#Function to check if a node is monophyletic 
checkMonophyly <- function(tree, node, mapping_df, tip_name, level) {
  offspring_leaves <- getOffspringLabels(tree, node)
  
  # Ensure that level exists in mapping_df
  if (!(level %in% colnames(mapping_df))) {
    stop(paste("The column", level, "does not exist in mapping_df"))
  }
  
  #Use mapping_df to determine groups
  offspring_groups <- mapping_df %>%
    filter(.data[[tip_name]] %in% offspring_leaves$label) %>%
    select(all_of(level)) %>%
    distinct()
  
  return(nrow(offspring_groups) == 1) #Returns TRUE if monophyletic
}

#Clade collapsing function (keeps tip with branch length closest to the median of
#corresponding monophyletic clade) 
collapse_clade <- function(tree, tips) {
  if (length(tips) < 2) {
    return(tree)  #Return unchanged if there’s only one tip
  }
  
  #Get branch lengths for the tips of provided clade
  tip_branches <- data.frame(
    tip = tree$tip.label,
    length = tree$edge.length[match(1:length(tree$tip.label), tree$edge[,2])]
  ) %>%
    filter(tip %in% tips)
  
  #Identify the tip with the longest branch
  #tip_to_keep <- tip_branches %>% 
  #  arrange(desc(length)) %>% 
  #  slice(1) %>% 
  #  pull(tip)
  
  #Identify the tip which length is closest to the median distance of the 
  #monophyletic clade
  
  tip_to_keep <- tip_branches %>% 
    mutate(distance_to_median = abs(length - median(length, na.rm = TRUE))) %>%
    arrange(distance_to_median) %>%
    slice(1) %>%
    pull(tip)
  
  #Remove all other tips except the longest branch
  tips_to_remove <- setdiff(tips, tip_to_keep)
  
  if (length(tips_to_remove) == 0) {
    message("No tips to drop for clade: ", paste(tips, collapse = ", "))
    return(tree)
  }
  
  tree <- drop.tip(tree, tips_to_remove)  # Drop tips except the longest one
  return(tree)
}

#Collapse function for the entire tree 
collapse_tree_to_clades <- function(tree, mapping_df, tip_name, level) {
  collapsed_tree <- tree
  
  # Compute node depths
  #node_depths <- node.depth.edgelength(tree)
  
  #Store internal nodes to iterate over
  internal_nodes <- (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode) #Total_Nodes_of_tree = length(tree$tip.label) + tree$Nnode
  #sorted_nodes <- internal_nodes[order(node_depths[internal_nodes])]  # Sort by depth (root first)
  
  for (node in internal_nodes) {
    
    if (is.na(node)) {
      warning("Node is NA. Skipping...")
      next  # Skip if the node is NA
    }
    
    #If offspring of node form monophyletic clade
    if (checkMonophyly(tree, node, mapping_df, tip_name, level)) {
      #Extract the leaves
      offspring_leaves <- getOffspringLabels(tree, node)$label
      #And collapse clade by only keeping the longest tip
      collapsed_tree <- collapse_clade(collapsed_tree, offspring_leaves)
    }
  }
  
  return(collapsed_tree)
}

#Rename tips of collapsed tree with taxonomy
rename_tree_tips <- function(collapsed_tree, mapping_df, tip_name, level) {
  #Rename tip labels based on the selected level of taxonomy
  new_tip_label <- mapping_df[match(collapsed_tree$tip.label, mapping_df[, tip_name]), level]
  
  #new_tip_label <- mapping_df %>%
  #  select(all_of(tip_name), all_of(level)) %>%
  #  right_join(tibble(tip = collapsed_tree$tip.label), by = setNames("tip", tip_name)) %>%
  #  pull(all_of(level))
  
  # Check for unmatched tips
  if (any(is.na(new_tip_label))) {
    missing_tips <- tree$tip.label[is.na(new_tip_label)]
    warning(paste0("The following tip(s) didn't match any entry in the reference dataframe: ",
                   paste(missing_tips, collapse = ", ")))
  }
  
  # Update the tree with the new tip labels
  collapsed_tree$tip.label <- new_tip_label
  return(collapsed_tree)
}

#Function to assign colors to tip and parent edges
create_colour_vector <- function(tree, mapping, tip_col, level_color) {
  edge_colors <- rep("black", nrow(tree$edge))
  
  # Get tip indices and handle missing values
  tip_indices <- match(tree$tip.label, mapping[[tip_col]])
  
  # Check for missing tips
  missing_tips <- setdiff(tree$tip.label, mapping[[tip_col]])
  if (length(missing_tips) > 0) {
    warning("Missing tips:", paste(missing_tips, collapse = ", "))
  }
  
  # Assign tip colors
  edge_tip_indices <- which(tree$edge[, 2] <= length(tree$tip.label))
  edge_colors[edge_tip_indices] <- mapping[[level_color]][tip_indices[tree$edge[edge_tip_indices, 2]]]
  
  # Propagate colors bottom-up
  for (i in nrow(tree$edge):1) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]
    
    if (edge_colors[i] == "black") {
      child_edges <- which(tree$edge[, 1] == child)
      child_colors <- unique(edge_colors[child_edges])
      
      if (length(child_colors) == 1 && child_colors != "black") {
        edge_colors[i] <- child_colors
      }
    }
  }
  
  return(edge_colors)
}
