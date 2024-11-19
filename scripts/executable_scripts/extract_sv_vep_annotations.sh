#!/usr/bin/env bash

set -o errexit -o nounset -o pipefail

vep_in=$1
cache_dir=$2
tmp=$3
vep_out=$4

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

python3 scripts/executable_scripts/extract_sv_vep_annotations.py ${tmp} ${vep_out}
