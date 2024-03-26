#!/usr/bin/env bash
# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
mamba activate rare_variants
python3.10 scripts/executable_scripts/train_test_predict_split_annotation.py \
--training /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/new-protein-lincRNA-100k-rare-position-aware-2.0/combined_annotation_pre_merge_gene_level_noimpute.tissue.tmm_corrected.no_vep.no_muscle.z3.csv \
--testings /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-1.0/combined_annotation_pre_merge_custom_maf_gene_level.no_control_gene.tmm_corrected.no_vep.csv \
--mode predict \
--min-af-impute-mode infer \
--out-prefix /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-1.0/tissue_specific_annotations_sv_pred_gene_level

python3.10 scripts/executable_scripts/train_test_predict_split_annotation.py \
--training /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/new-protein-lincRNA-100k-rare-position-aware-2.0/combined_annotation_pre_merge_gene_level_noimpute.tissue.tmm_corrected.no_vep.no_muscle.z3.csv \
--testings /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-1.0/combined_annotation_pre_merge_custom_maf.no_control_gene.tmm_corrected.no_vep.csv \
--mode predict \
--min-af-impute-mode infer --out-prefix \
/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN-protein-lincRNA-100k-VNTR-rare-position-aware-1.0/tissue_specific_annotations_sv_pred
