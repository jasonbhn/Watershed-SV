#!/usr/bin/env bash
# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
mamba activate rare_variants
cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/repeat_expansion_testing/UDN.VNTR.mean_neighbor_distance.kneighbors_10.extreme_outliers.RARE.vcf /tmp/
cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/finalized_UDN_expression/combined_Blood_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.noGlobalOutliers.zscore.tsv /tmp/
cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/finalized_UDN_expression/combined_Fibroblast_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.noGlobalOutliers.zscore.tsv /tmp/
cp -r UDN-protein-lincRNA-10k-VNTR-final-version /tmp/ 
python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/UDN.VNTR.mean_neighbor_distance.kneighbors_10.extreme_outliers.RARE.vcf \
--genotypes /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/genes.bed \
--gene-sv /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/gene_sv.10000.bed \
--annotation-dir /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/ \
--outfile UDN-protein-lincRNA-10k-VNTR-final-version/combined_annotation_pre_merge_custom_maf.no_control_gene.tmm_corrected.Blood.csv \
--expressions /tmp/combined_Blood_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.noGlobalOutliers.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/custom_UDN_maf.tsv \
--length-mode upload-VNTR \
--length-file /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/custom_UDN_length.tsv \
--CN-mode extract \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 10000 \

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/UDN.VNTR.mean_neighbor_distance.kneighbors_10.extreme_outliers.RARE.vcf \
--genotypes /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/genes.bed \
--gene-sv /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/gene_sv.10000.bed \
--annotation-dir /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/ \
--outfile UDN-protein-lincRNA-10k-VNTR-final-version/combined_annotation_pre_merge_custom_maf.no_control_gene.tmm_corrected.Fibroblast.csv \
--expressions /tmp/combined_Fibroblast_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.noGlobalOutliers.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/custom_UDN_maf.tsv \
--length-mode upload-VNTR \
--length-file /tmp/UDN-protein-lincRNA-10k-VNTR-final-version/intermediates/custom_UDN_length.tsv \
--CN-mode extract \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 10000 \

