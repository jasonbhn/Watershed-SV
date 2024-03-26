#!/bin/bash
vep_path=/home/bni1/vep.sif 
vep_in=$1
cache_dir=/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/vep_cache 
tmp=/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN022194_Clair3_vcf/tmp.csv 

# sort the input vep. 
 vep \
-i ${vep_in} \
--dir_cache ${cache_dir} \
-o ${tmp} \
--format vcf \
--tab \
--fork 4 \
--cache \
--force \
--regulatory \
--variant_class \
--max_af \
--af_gnomadg

