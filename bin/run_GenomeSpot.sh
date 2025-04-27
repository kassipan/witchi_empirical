#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate GenomeSPOT

file=$1

cd /local/two/panag007/software/GenomeSPOT
python -m genome_spot.genome_spot --models models --contigs $file --proteins $(echo $file | sed 's/fna/faa/g') --output predictions/$(basename $file .fna)
