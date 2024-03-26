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
#SBATCH --time=1:40:0
#SBATCH --account=smontgom
#SBATCH --mem=200G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
module load anaconda
eval "$(conda shell.bash hook)"
conda activate rare_variants

python3.9 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf UDN_dataset/UDN.cohort_combined.GRCh38.consensus_SVs.jasmine_irisRefined.major_chr_seq_resolved.paragraph_1kgp_maf.reheader.annotated.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-10.0/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-10.0/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-10.0/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-10.0/intermediates/ \
--outfile UDN-protein-lincRNA-10k-10.0/combined_annotation_Blood_custom_maf_new_outlier.csv \
--expressions UDN_dataset/UDN_Blood_regress_top_60_PCs_Sex_no_global_outlier.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file UDN-protein-lincRNA-10k-10.0/intermediates/custom_UDN_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000

python3.9 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf UDN_dataset/UDN.cohort_combined.GRCh38.consensus_SVs.jasmine_irisRefined.major_chr_seq_resolved.paragraph_1kgp_maf.reheader.annotated.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-10.0/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-10.0/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-10.0/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-10.0/intermediates/ \
--outfile UDN-protein-lincRNA-10k-10.0/combined_annotation_Fibroblast_custom_maf_new_outlier.csv \
--expressions UDN_dataset/UDN_Fibroblast_regress_top_60_PCs_Sex_no_global_outlier.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file UDN-protein-lincRNA-10k-10.0/intermediates/custom_UDN_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000
