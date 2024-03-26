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
#SBATCH --time=4:0:0
#SBATCH --account=smontgom
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
eval "$(conda shell.bash hook)"
conda activate rare_variants

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf input/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes new-protein-lincRNA-10k-9.0/intermediates/pipeline_input_genotypes.tsv \
--genes new-protein-lincRNA-10k-9.0/intermediates/genes.bed \
--gene-sv new-protein-lincRNA-10k-9.0/intermediates/gene_sv.10000.bed \
--annotation-dir new-protein-lincRNA-10k-9.0/intermediates/ \
--outfile new-protein-lincRNA-10k-9.0/combined_annotation_pre_merge_gene_level_noimpute_blood.csv \
--expressions UDN_dataset/finalized_UDN_expression/combined_Blood_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file new-protein-lincRNA-10k-9.0/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file new-protein-lincRNA-10k-9.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 10000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf input/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes new-protein-lincRNA-10k-9.0/intermediates/pipeline_input_genotypes.tsv \
--genes new-protein-lincRNA-10k-9.0/intermediates/genes.bed \
--gene-sv new-protein-lincRNA-10k-9.0/intermediates/gene_sv.10000.bed \
--annotation-dir new-protein-lincRNA-10k-9.0/intermediates/ \
--outfile new-protein-lincRNA-10k-9.0/combined_annotation_pre_merge_noimpute_blood.csv \
--expressions UDN_dataset/finalized_UDN_expression/combined_Blood_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file new-protein-lincRNA-10k-9.0/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file new-protein-lincRNA-10k-9.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 10000
