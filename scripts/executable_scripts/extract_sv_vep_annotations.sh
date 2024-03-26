#!/bin/bash
vep_path=$1
vep_in=$2
cache_dir=$3
tmp=$4
vep_out=$5
flank=$6
# sort the input vep. 
module load ensembl-vep/109.3

sort -k1,1 -k2,2n ${vep_in} | vep \
-o ${tmp} \
--format ensembl \
--verbose \
--cache \
--dir ${cache_dir} \
--tab \
--fields "Uploaded_variation,Gene,Feature_type,Consequence,IMPACT" \
--fork 4 \
--force \
--regulatory \
--overlaps \
--distance 10000
module unload ensembl-vep/109.3
python3.10 scripts/executable_scripts/extract_sv_vep_annotations.py ${tmp} ${vep_out}
