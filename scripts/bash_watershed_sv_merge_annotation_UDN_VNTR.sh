#!/usr/bin/env bash
# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
mamba activate rare_variants
cp  /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/repeat_expansion_testing/UDN.VNTR.mean_neighbor_distance.kneighbors_10.extreme_outliers.RARE.vcf /tmp/
cp  /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/input/UDN_LRS_ABC_udn_expression.tmm_corrected.csv /tmp/
cp -r UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0 /tmp/ 
python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/UDN.VNTR.mean_neighbor_distance.kneighbors_10.extreme_outliers.RARE.vcf \
--genotypes /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/genes.bed \
--gene-sv /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/ \
--outfile UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/combined_annotation_pre_merge_custom_maf.no_control_gene.tmm_corrected.csv \
--expressions /tmp/UDN_LRS_ABC_udn_expression.tmm_corrected.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/custom_UDN_maf.tsv \
--length-mode upload-VNTR \
--length-file /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/custom_UDN_length.tsv \
--CN-mode extract \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000 \
--minimum-support-tissue-count 1

python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/UDN.VNTR.mean_neighbor_distance.kneighbors_10.extreme_outliers.RARE.vcf  \
--genotypes /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/genes.bed \
--gene-sv /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/ \
--outfile UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/combined_annotation_pre_merge_custom_maf_gene_level.no_control_gene.tmm_corrected.csv \
--expressions /tmp/UDN_LRS_ABC_udn_expression.tmm_corrected.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/custom_UDN_maf.tsv \
--length-mode upload-VNTR \
--length-file /tmp/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-2.0/intermediates/custom_UDN_length.tsv \
--CN-mode extract \
--collapse-mode gene \
--remove-control-genes \
--flank 100000 \
--minimum-support-tissue-count 1
