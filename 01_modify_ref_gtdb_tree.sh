#!/bin/bash

REF_TREE=data/ar53_r220.tree.renamed
WITCHI_TREE=data/ar53_msa_reps_r220_squared_s10_pruned60.fasta.PMSF_LGC10G4F.main.treefile.renamed.noSHaLRT
mkdir -p trees

#Remove node labels from ref tree
sed -E "s/'([^:]+):[^']+'/\1/g" $REF_TREE > trees/ar53_r220.tree.renamed.mod

nw_reroot -s trees/ar53_r220.tree.renamed.mod $(nw_labels -I trees/ar53_r220.tree.renamed.mod | grep "p__Altiarchaeota\|p__EX4484-52") \
> trees/ar53_r220.tree.renamed.mod.dpann.root

nw_reroot -s $WITCHI_TREE $(nw_labels -I $WITCHI_TREE | grep "p__Altiarchaeota\|p__EX4484-52") \
> trees/ar53_msa_reps_r220_squared_s10_pruned60.fasta.PMSF_LGC10G4F.main.treefile.renamed.noSHaLRT.dpann.root

#Extract p__Halobacteriota clade
nw_clade -r trees/ar53_r220.tree.renamed.mod '.*p__Halobacteriota.*' > trees/ar53_r220_p__Halobacteriota_subtree.tree
nw_clade -r $WITCHI_TREE '.*p__Halobacteriota.*' > trees/ar53_r220_squared_s10_pruned60_p__Halobacteriota_subtree.tree


