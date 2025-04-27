#!/usr/bin/env Rscript

# Load libraries & functions
library(ape)
library(phytools)
library(dplyr)
library(tidyr)
library(treeio)
library(tidytree)
print("sourcing functions.R")
source("bin/functions.R")
print("done sourcing")

full_tree_path <- "trees/ar53_r220.tree.renamed.mod"
pruned_tree_path <- 
  "trees/ar53_msa_reps_r220_squared_s10_pruned60.fasta.PMSF_LGC10G4F.main.treefile.renamed.noSHaLRT.dpann.root"
collapse_tsv <- "data/gtdb_r220.map"

#Load trees
tree_full <- read.tree(full_tree_path)
tree_full <- midpoint.root(tree_full)
tree_pruned <- read.tree(pruned_tree_path)

#Load collapse mapping file
collapse_df <- read.table(collapse_tsv, sep = "\t", header = TRUE, 
                          stringsAsFactors = FALSE)

################### Class-level Tanglegram (supergroup color) ###################
#Collapse both trees based on desired level of taxonomy
tree_full_collapsed <- collapse_tree_to_clades(tree = tree_full, mapping_df = collapse_df, tip_name = "tip_label", level = "class")
tree_full_collapsed_renamed <- rename_tree_tips(tree_full_collapsed, collapse_df, tip_name = "tip_label", level = "class")
tree_full_collapsed_renamed$tip.label <- gsub("^c__", "", tree_full_collapsed_renamed$tip.label)

tree_pruned_collapsed <- collapse_tree_to_clades(tree = tree_pruned, mapping_df = collapse_df, tip_name = "tip_label", level = "class")
tree_pruned_collapsed_renamed <- rename_tree_tips(tree_pruned_collapsed, collapse_df, tip_name = "tip_label", level = "class")
tree_pruned_collapsed_renamed$tip.label <- gsub("^c__", "", tree_pruned_collapsed_renamed$tip.label)

#plot.phylo(tree_full_collapsed_renamed, type = "phylogram", edge.width = 1, show.node.label = TRUE)
#title("Collapsed Full Tree")
#add.scale.bar(length = 0.1) 

#plot.phylo(tree_pruned_collapsed_renamed, type = "phylogram", edge.width = 1, show.node.label = TRUE)
#title("Collapsed Pruned Tree")

#Association matrix for the two trees
association <- cbind(tree_full_collapsed_renamed$tip.label, tree_full_collapsed_renamed$tip.label)

#Draw tanglegram
trees_cophylo <- cophylo(tree_full_collapsed_renamed, tree_pruned_collapsed_renamed, 
                         assoc = association, rotate = TRUE, print = TRUE)

collapse_df$class <- gsub("^c__", "", collapse_df$class)

#Assign colors to tip and parent edges in both left and right trees of cophylo
edge_colour_left <- create_colour_vector(trees_cophylo$trees[[1]], collapse_df, "class", "supergroup_color")
edge_colour_right <- create_colour_vector(trees_cophylo$trees[[2]], collapse_df, "class", "supergroup_color")
edge.col <- list(left = edge_colour_left, right = edge_colour_right)
link.col <- as.character(collapse_df$supergroup_color[match(trees_cophylo$assoc[,1], collapse_df$class)])


pdf("figures/tanglegram_arGTDB_vs_witchi_class_sepergroupColor.pdf", width = 9, height = 11)
plot(trees_cophylo,
     link.type = "curved",
     edge.col = edge.col,
     fsize = 1,
     link.lwd = 3,
     link.lty = "solid",
     link.col = link.col,
     tip.len = 0,
     lwd = 2.3,
     pts=FALSE,
     scale.bar=rep(0.1,2))

#Extract bootstrap values from both trees
bootstrap_left <- as.numeric(trees_cophylo$trees[[1]]$node.label)
bootstrap_right <- as.numeric(trees_cophylo$trees[[2]]$node.label)

#Set bootstrap values to black if Ufboot >95, else NA
bootstrap_left[bootstrap_left > 95] <- "black"
bootstrap_left[bootstrap_left < 95] <- NA
bootstrap_right[bootstrap_right > 95] <- "black"
bootstrap_right[bootstrap_right < 95] <- NA

valid_nodes_left <- which(!is.na(bootstrap_left) & bootstrap_left > 95) + length(trees_cophylo$trees[[1]]$tip.label)
valid_nodes_right <- which(!is.na(bootstrap_right) & bootstrap_right > 95) + length(trees_cophylo$trees[[2]]$tip.label)

# Add bootstrap values as node labels only for valid nodes
nodelabels.cophylo(node = valid_nodes_left,
                   pch = 21,  
                   frame = "circle",
                   which = "left",
                   cex = 1.1,  
                   bg = "black",  
                   col = "black")

nodelabels.cophylo(node = valid_nodes_right,
                   pch = 21,  
                   frame = "circle",
                   which = "right",
                   cex = 1.1,  
                   bg = "black",  
                   col = "black" )

legend_data <- unique(collapse_df[, c("supergroup", "supergroup_color")])

legend("bottom", 
       legend = legend_data$supergroup, 
       col = legend_data$supergroup_color, 
       bty = "n", 
       lty = 1, 
       cex = 0.8, 
       lwd = 5, 
       ncol = 2,  # Set to two columns
       text.font = 1,
       xpd = TRUE)

legend("bottomleft",
       legend = "≥95% UFBoot Support",
       pch = 21,  
       pt.bg = "black",  
       col = "black",  
       pt.cex = 1.1,  # Same as nodelabels.cophylo cex
       bty = "n")

dev.off()

################### Class-level Tanglegram (phylum color) ###################

#Assign colors to tip and parent edges in both left and right trees of cophylo
edge_colour_left <- create_colour_vector(trees_cophylo$trees[[1]], collapse_df, "class", "phylum_color")
edge_colour_right <- create_colour_vector(trees_cophylo$trees[[2]], collapse_df, "class", "phylum_color")
edge.col <- list(left = edge_colour_left, right = edge_colour_right)
link.col <- as.character(collapse_df$phylum_color[match(trees_cophylo$assoc[,1], collapse_df$class)])

pdf("figures/tanglegram_arGTDB_vs_witchi_class_phylumColor.pdf", width = 8.5, height = 11)
plot(trees_cophylo,
     link.type = "curved",
     edge.col = edge.col,
     fsize = 1,
     link.lwd = 4,
     link.lty = "solid",
     link.col = link.col,
     tip.len = 0,
     lwd = 3,
     pts=FALSE)

#Extract bootstrap values from both trees
bootstrap_left <- as.numeric(trees_cophylo$trees[[1]]$node.label)
bootstrap_right <- as.numeric(trees_cophylo$trees[[2]]$node.label)

#Set bootstrap values to black if Ufboot >95, else NA
bootstrap_left[bootstrap_left > 95] <- "black"
bootstrap_left[bootstrap_left < 95] <- NA
bootstrap_right[bootstrap_right > 95] <- "black"
bootstrap_right[bootstrap_right < 95] <- NA
    
valid_nodes_left <- which(!is.na(bootstrap_left) & bootstrap_left > 95) + length(trees_cophylo$trees[[1]]$tip.label)
valid_nodes_right <- which(!is.na(bootstrap_right) & bootstrap_right > 95) + length(trees_cophylo$trees[[2]]$tip.label)
    
#Add bootstrap values as node labels only for valid nodes
nodelabels.cophylo(node = valid_nodes_left,
                       pch = 21,  
                       frame = "circle",
                       which = "left",
                       cex = 1,  
                       bg = "black",  
                       col = "black")
    
nodelabels.cophylo(node = valid_nodes_right,
                       pch = 21,  
                       frame = "circle",
                       which = "right",
                       cex = 1,  
                       bg = "black",  
                       col = "black" )

legend_data <- unique(collapse_df[, c("phylum", "phylum_color")])

legend("bottom", 
       legend = legend_data$phylum, 
       col = legend_data$phylum_color, 
       bty = "n", 
       lty = 1, 
       cex = 0.8, 
       lwd = 5, 
       ncol = 2,  # Set to two columns
       text.font = 1,
       xpd = TRUE)  

dev.off()
    
########################### Phylum-level Tanglegram ############################
tree_full_collapsed <- collapse_tree_to_clades(tree = tree_full, mapping_df = collapse_df, tip_name = "tip_label", level = "phylum")
tree_full_collapsed_renamed <- rename_tree_tips(tree_full_collapsed, collapse_df, tip_name = "tip_label", level = "phylum")
tree_full_collapsed_renamed$tip.label <- gsub("^p__", "", tree_full_collapsed_renamed$tip.label)

tree_pruned_collapsed <- collapse_tree_to_clades(tree = tree_pruned, mapping_df = collapse_df, tip_name = "tip_label", level = "phylum")
tree_pruned_collapsed_renamed <- rename_tree_tips(tree_pruned_collapsed, collapse_df, tip_name = "tip_label", level = "phylum")
tree_pruned_collapsed_renamed$tip.label <- gsub("^p__", "", tree_pruned_collapsed_renamed$tip.label)

#plot.phylo(tree_full_collapsed_renamed, type = "phylogram", edge.width = 1, show.node.label = TRUE)
#title("Collapsed Full Tree")
#add.scale.bar(length = 0.1) 
    
#plot.phylo(tree_pruned_collapsed_renamed, type = "phylogram", edge.width = 1, show.node.label = TRUE)
#title("Collapsed Pruned Tree")
    
#Association matrix for the two trees
association <- cbind(tree_full_collapsed_renamed$tip.label, tree_full_collapsed_renamed$tip.label)
    
#Draw tanglegram
trees_cophylo <- cophylo(tree_full_collapsed_renamed, tree_pruned_collapsed_renamed, 
                             assoc = association, rotate = TRUE, print = TRUE)
    
collapse_df$phylum <- gsub("^p__", "", collapse_df$phylum)
    
#Assign colors to tip and parent edges in both left and right trees of cophylo
edge_colour_left <- create_colour_vector(trees_cophylo$trees[[1]], collapse_df, "phylum", "supergroup_color")
edge_colour_right <- create_colour_vector(trees_cophylo$trees[[2]], collapse_df, "phylum", "supergroup_color")
edge.col <- list(left = edge_colour_left, right = edge_colour_right)
link.col <- as.character(collapse_df$supergroup_color[match(trees_cophylo$assoc[,1], collapse_df$phylum)])

    
    
pdf("figures/tanglegram_arGTDB_vs_witchi_phylum_supergroupColor.pdf", width = 9, height = 7)
plot(trees_cophylo,
     link.type = "curved",
     edge.col = edge.col,
     fsize = 1,
     link.lwd = 2.5,
     link.lty = "solid",
     link.col = link.col,
     tip.len = 0,
     lwd = 2.8,
     pts=FALSE,
     scale.bar=rep(0.1,2))
    
#Extract bootstrap values from both trees
bootstrap_left <- as.numeric(trees_cophylo$trees[[1]]$node.label)
bootstrap_right <- as.numeric(trees_cophylo$trees[[2]]$node.label)
    
#Set bootstrap values to black if Ufboot >95, else NA
bootstrap_left[bootstrap_left > 95] <- "black"
bootstrap_left[bootstrap_left < 95] <- NA
bootstrap_right[bootstrap_right > 95] <- "black"
bootstrap_right[bootstrap_right < 95] <- NA
        
valid_nodes_left <- which(!is.na(bootstrap_left) & bootstrap_left > 95) + length(trees_cophylo$trees[[1]]$tip.label)
valid_nodes_right <- which(!is.na(bootstrap_right) & bootstrap_right > 95) + length(trees_cophylo$trees[[2]]$tip.label)
        
# Add bootstrap values as node labels only for valid nodes
nodelabels.cophylo(node = valid_nodes_left,
                   pch = 21,
                   frame = "circle",
                   which = "left",
                   cex = 1.3,  
                   bg = "black",
                   col = "black")
        
nodelabels.cophylo(node = valid_nodes_right,
                   pch = 21,  
                   frame = "circle",
                   which = "right",
                   cex = 1.3,  
                   bg = "black",  
                   col = "black" )
        
legend("bottomleft",
       legend = "≥95% UFBoot Support",
       pch = 21,  
       pt.bg = "black",  
       col = "black",  
       pt.cex = 1.3,  #Same as nodelabels.cophylo cex
       bty = "n")

legend_data <- unique(collapse_df[, c("supergroup", "supergroup_color")])

legend("bottom", 
       legend = legend_data$supergroup, 
       col = legend_data$supergroup_color, 
       bty = "n", 
       lty = 1, 
       cex = 0.8, 
       lwd = 5, 
       ncol = 2,
       text.font = 1,
       xpd = TRUE)  

dev.off()

################### Halobacteriota Tanglegram (phylum color) ###################

halo_ref_path <- "trees/ar53_r220_p__Halobacteriota_subtree.tree"
halo_pruned_path <- "trees/ar53_r220_squared_s10_pruned60_p__Halobacteriota_subtree.tree"

halo_ref <- read.tree(halo_ref_path)
halo_pruned <- read.tree(halo_pruned_path)

halo_ref_collapsed <- collapse_tree_to_clades(tree = halo_ref, mapping_df = collapse_df, tip_name = "tip_label", level = "class")
halo_ref_collapsed_renamed <- rename_tree_tips(halo_ref_collapsed, collapse_df, tip_name = "tip_label", level = "class")
halo_ref_collapsed_renamed$tip.label <- gsub("^c__", "", halo_ref_collapsed_renamed$tip.label)

halo_pruned_collapsed <- collapse_tree_to_clades(tree = halo_pruned, mapping_df = collapse_df, tip_name = "tip_label", level = "class")
halo_pruned_collapsed_renamed <- rename_tree_tips(halo_pruned_collapsed, collapse_df, tip_name = "tip_label", level = "class")
halo_pruned_collapsed_renamed$tip.label <- gsub("^c__", "", halo_pruned_collapsed_renamed$tip.label)

#plot.phylo(halo_ref_collapsed_renamed, type = "phylogram", edge.width = 1, show.node.label = TRUE)
#title("Collapsed Halo Tree Ref")
#add.scale.bar(length = 0.1) 

#plot.phylo(halo_pruned_collapsed_renamed, type = "phylogram", edge.width = 1, show.node.label = TRUE)
#title("Collapsed Halo Tree Pruned")

#Association matrix for the two trees
association_halo <- cbind(halo_ref_collapsed_renamed$tip.label, halo_ref_collapsed_renamed$tip.label)

#Draw tanglegram
halo_trees_cophylo <- cophylo(halo_ref_collapsed_renamed, halo_pruned_collapsed_renamed, 
                         assoc = association_halo, rotate = TRUE, print = TRUE)

collapse_df$class <- gsub("^c__", "", collapse_df$class)

#Assign colors to tip and parent edges in both left and right trees of cophylo
edge_colour_left <- create_colour_vector(halo_trees_cophylo$trees[[1]], collapse_df, "class", "supergroup_color")
edge_colour_right <- create_colour_vector(halo_trees_cophylo$trees[[2]], collapse_df, "class", "supergroup_color")
edge.col <- list(left = edge_colour_left, right = edge_colour_right)
link.col <- as.character(collapse_df$supergroup_color[match(halo_trees_cophylo$assoc[,1], collapse_df$class)])

pdf("figures/tanglegram_arGTDB_vs_witchi_Halobacteriota.pdf", width = 9, height = 5)
plot(halo_trees_cophylo,
     link.type = "curved",
     edge.col = edge.col,
     fsize = 1,
     link.lwd = 3,
     link.lty = "solid",
     link.col = link.col,
     tip.len = 0,
     lwd = 2.8,
     pts=FALSE,
     scale.bar=rep(0.1,2))

#Extract bootstrap values from both trees
bootstrap_left <- as.numeric(halo_trees_cophylo$trees[[1]]$node.label)
bootstrap_right <- as.numeric(halo_trees_cophylo$trees[[2]]$node.label)

#Set bootstrap values to black if Ufboot >95, else NA
bootstrap_left[bootstrap_left > 95] <- "black"
bootstrap_left[bootstrap_left < 95] <- NA
bootstrap_right[bootstrap_right > 95] <- "black"
bootstrap_right[bootstrap_right < 95] <- NA
    
valid_nodes_left <- which(!is.na(bootstrap_left) & bootstrap_left > 95) + length(halo_trees_cophylo$trees[[1]]$tip.label)
valid_nodes_right <- which(!is.na(bootstrap_right) & bootstrap_right > 95) + length(halo_trees_cophylo$trees[[2]]$tip.label)
    
# Add bootstrap values as node labels only for valid nodes
nodelabels.cophylo(node = valid_nodes_left,
                   pch = 21,
                   frame = "circle",
                   which = "left",
                   cex = 1.3,
                   bg = "black",
                   col = "black")
    
nodelabels.cophylo(node = valid_nodes_right,
                   pch = 21,
                   frame = "circle",
                   which = "right",
                   cex = 1.3,
                   bg = "black",
                   col = "black" )
    
legend("bottomleft",
       legend = "≥95% UFBoot Support",
       pch = 21,  
       pt.bg = "black",  
       col = "black",  
       pt.cex = 1.3,  # Same as nodelabels.cophylo cex
       bty = "n")

dev.off()

    
