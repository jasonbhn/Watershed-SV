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
cp -r new-protein-lincRNA-100k-rare-baseline /tmp/
echo "here"
python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_gene_level_noimpute.tissue.tmm_corrected.z3.csv \
--expressions /tmp/UDN_LRS_ABC_gtex_expression.tmm_corrected.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 100000 \

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_noimpute.tissue.tmm_corrected.z3.csv \
--expressions /tmp/UDN_LRS_ABC_gtex_expression.tmm_corrected.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000 \
