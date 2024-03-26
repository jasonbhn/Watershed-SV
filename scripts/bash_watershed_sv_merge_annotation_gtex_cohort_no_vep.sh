#!/usr/bin/env bash
# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
mamba activate rare_variants
cp input/gtex.lumpy.gs.melt.high_conf.vcf.gz /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz
#cp input/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.z3.txt /tmp/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.z3.txt
cp input/UDN_LRS_ABC_gtex_expression.tmm_corrected.csv /tmp/UDN_LRS_ABC_gtex_expression.tmm_corrected.csv
cp -r new-protein-lincRNA-100k-rare-position-aware-2.0 /tmp/
echo "here"
python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/ \
--outfile new-protein-lincRNA-100k-rare-position-aware-2.0/combined_annotation_pre_merge_gene_level_noimpute.tissue.tmm_corrected.no_vep.z3.csv \
--expressions /tmp/UDN_LRS_ABC_gtex_expression.tmm_corrected.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 100000 \
--zscore-threshold 3

python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/ \
--outfile new-protein-lincRNA-100k-rare-position-aware-2.0/combined_annotation_pre_merge_noimpute.tissue.tmm_corrected.no_vep.z3.csv \
--expressions /tmp/UDN_LRS_ABC_gtex_expression.tmm_corrected.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000 \
--zscore-threshold 3
