#!/usr/bin/env bash
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
#SBATCH --time=3:0:0
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --account=smontgom
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
mamba activate rare_variants

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf UDN_dataset/UDN.VNTR.rare_extreme_only.Watershed.vcf \
--genotypes UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/ \
--outfile UDN-protein-lincRNA-10k-VNTR-2.0/combined_annotation_PCA_expression_Blood_custom_maf.csv \
--expressions UDN_dataset/finalized_UDN_expression/combined_Blood_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/custom_UDN_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000
python3.10 scripts/executable_scripts/train_test_predict_split_annotation.py --training new-protein-lincRNA-10k-9.0/combined_annotation_pre_merge_gene_level_noimpute.csv --mode predict --out-prefix UDN-protein-lincRNA-10k-VNTR-2.0/Blood_annotations_custom_maf_PCA_expression --min-af-impute-mode infer --testings UDN-protein-lincRNA-10k-VNTR-2.0/combined_annotation_PCA_expression_Blood_custom_maf.csv
python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf UDN_dataset/UDN.VNTR.rare_extreme_only.Watershed.vcf \
--genotypes UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/ \
--outfile UDN-protein-lincRNA-10k-VNTR-2.0/combined_annotation_PCA_expression_Fibroblast_custom_maf.csv \
--expressions UDN_dataset/finalized_UDN_expression/combined_Fibroblast_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file UDN-protein-lincRNA-10k-VNTR-2.0/intermediates/custom_UDN_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000
python3.10 scripts/executable_scripts/train_test_predict_split_annotation.py --training new-protein-lincRNA-10k-9.0/combined_annotation_pre_merge_gene_level_noimpute.csv --mode predict --out-prefix UDN-protein-lincRNA-10k-VNTR-2.0/Fibroblast_annotations_custom_maf_PCA_expression --min-af-impute-mode infer --testings UDN-protein-lincRNA-10k-VNTR-2.0/combined_annotation_PCA_expression_Fibroblast_custom_maf.csv
