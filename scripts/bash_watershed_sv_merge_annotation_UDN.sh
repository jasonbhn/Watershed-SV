#!/usr/bin/env bash
# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
mamba activate rare_variants
cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN_MergedSV_vcfs/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.noBNDTRA.vcf.gz /tmp/ 
cp input/UDN_LRS_ABC_udn_expression.tmm_corrected.csv /tmp/
cp -r UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/ /tmp/ 
echo "here"
python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.noBNDTRA.vcf.gz \
--genotypes /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/intermediates/genes.bed \
--gene-sv /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/ \
--outfile UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/combined_annotation_pre_merge_custom_maf.no_control_gene.tmm_corrected.no_vep.csv \
--expressions /tmp/UDN_LRS_ABC_udn_expression.tmm_corrected.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode extract \
--maf-field ALNG_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000 \
--minimum-support-tissue-count 1


python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.noBNDTRA.vcf.gz \
--genotypes /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/intermediates/genes.bed \
--gene-sv /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/ \
--outfile UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/combined_annotation_pre_merge_custom_maf_gene_level.no_control_gene.tmm_corrected.no_vep.csv \
--expressions /tmp/UDN_LRS_ABC_udn_expression.tmm_corrected.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode extract \
--maf-field ALNG_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene \
--remove-control-genes \
--flank 100000 \
--minimum-support-tissue-count 1


