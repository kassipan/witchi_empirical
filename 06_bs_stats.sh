#!/bin/bash

for tree in trees/*newick; do
    echo $tree
    phykit bipartition_support_stats $tree
done
