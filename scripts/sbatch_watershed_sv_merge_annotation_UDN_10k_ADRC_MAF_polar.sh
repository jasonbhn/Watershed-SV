#!/bin/bash
# ---------------------------------------------------
# The Advanced Research Computing at Hopkins (ARCH)
# User and Application Support < help@rockfish.jhu.edu >
#
# SLURM script to run the JupyterLab
#
# ---------------------------------------------------
#  INPUT ENVIRONMENT VARIABLES
# ---------------------------------------------------
#SBATCH --job-name=WtsdSV
#SBATCH --time=2:0:0
#SBATCH --account=smontgom
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
eval "$(conda shell.bash hook)"
conda activate rare_variants

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN_MergedSV_vcfs/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.rare.noBNDTRA.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-ADRC-new/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-ADRC-new/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-ADRC-new/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-ADRC-new/intermediates/ \
--outfile UDN-protein-lincRNA-10k-ADRC-new/combined_annotation_pre_merge_gene_Blood_custom_maf_gene_level.csv \
--expressions UDN_dataset/finalized_UDN_expression/combined_Blood_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode extract \
--maf-field ALNG_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene \
--flank 10000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN_MergedSV_vcfs/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.rare.noBNDTRA.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-ADRC-new/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-ADRC-new/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-ADRC-new/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-ADRC-new/intermediates/ \
--outfile UDN-protein-lincRNA-10k-ADRC-new/combined_annotation_pre_merge_gene_Blood_custom_maf.csv \
--expressions UDN_dataset/finalized_UDN_expression/combined_Blood_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode extract \
--maf-field ALNG_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN_MergedSV_vcfs/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.rare.noBNDTRA.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-ADRC-new/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-ADRC-new/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-ADRC-new/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-ADRC-new/intermediates/ \
--outfile UDN-protein-lincRNA-10k-ADRC-new/combined_annotation_pre_merge_gene_Fibroblast_custom_maf_gene_level.csv \
--expressions UDN_dataset/finalized_UDN_expression/combined_Fibroblast_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode extract \
--maf-field ALNG_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene \
--flank 10000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN_MergedSV_vcfs/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.rare.noBNDTRA.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-ADRC-new/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-ADRC-new/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-ADRC-new/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-ADRC-new/intermediates/ \
--outfile UDN-protein-lincRNA-10k-ADRC-new/combined_annotation_pre_merge_gene_Fibroblast_custom_maf.csv \
--expressions UDN_dataset/finalized_UDN_expression/combined_Fibroblast_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode extract \
--maf-field ALNG_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000
