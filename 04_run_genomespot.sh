#!/bin/bash

#Archaea species representatives of GTDB r220 as input genomes
parallel -j 20 "prodigal -i {} -a data/faa/{/.}.faa -p meta" ::: data/fna/*.fna

cd /local/two/panag007/software/GenomeSPOT

ls data/fna/*.fna | parallel -j20 "bash bin/run_GenomeSpot.sh"

cd /local/two/panag007/15.wichi/empirical_datasets/gtdb_r220/genomes_reps
mkdir -p data/genomespot_predictions
mv /local/two/panag007/software/GenomeSPOT/predictions/*tsv genomespot_predictions

#Make taxonomy mapping file
echo -e "filename\tncbi_accession" > data/accesions.map
ls faa/*faa | sed 's/.*\///g' | sed 's/.faa//g' | awk -F"\t" '{split($1, arr, "_"); print $1 "\t" arr[1] "_" arr[2]}' >> data/accesions.map

Rscript --vanilla /local/two/panag007/scripts/left_join_files.R data/accesions.map data/ar53_metadata_r220.reps.tsv.map ncbi_accession data/taxonomy_map.tsv
