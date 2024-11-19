#!/usr/bin/env bash

set -o nounset -o pipefail -o errexit

gene_sv=$1
split_roadmap=$2
output=$3
state=$4

bedtools intersect -wa -wb -a $gene_sv -b $split_roadmap/*.$state.bed | 
sort -k4,4 -k6,6 | 
bedtools groupby -i stdin -g 4,6 -c 7 -o count_distinct | 
awk -v s=$state 'BEGIN{printf "SV\tGene\tstate_%s\n",s};{print}' >$output

