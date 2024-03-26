#!/usr/bin/env bash
# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate rare_variants

cp -r new-protein-lincRNA-100k-rare-baseline /tmp/
cp collapse_annotation_instructions.tsv /tmp/
python3.10 scripts/executable_scripts/eval_watershed_prep.py \
--gene-sv-annotation /tmp/new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_noimpute.tissue.tmm_corrected.z3.csv \
--gene-annotation /tmp/new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_gene_level_noimpute.tissue.tmm_corrected.z3.csv \
--collapse-instructions /tmp/collapse_annotation_instructions.tsv \
--pval-threshold 0.0027 \
--output new-protein-lincRNA-100k-rare-baseline/eval_watarshed_ABC_regional.tissue_tmm.z3.tsv
