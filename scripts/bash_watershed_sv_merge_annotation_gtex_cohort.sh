#!/usr/bin/env bash
cp input/gtex.lumpy.gs.melt.high_conf.vcf.gz /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz
cp input/UDN_LRS_ABC_gtex_expression.csv /tmp/UDN_LRS_ABC_gtex_expression.csv
cp -r new-protein-lincRNA-100k-rare-position-aware /tmp/

python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-position-aware/ \
--outfile new-protein-lincRNA-100k-rare-position-aware/combined_annotation_pre_merge_gene_level_noimpute.expression.csv \
--expressions /tmp/UDN_LRS_ABC_gtex_expression.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 100000

python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-position-aware/ \
--outfile new-protein-lincRNA-100k-rare-position-aware/combined_annotation_pre_merge_noimpute.expression.csv \
--expressions /tmp/UDN_LRS_ABC_gtex_expression.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-position-aware/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000

