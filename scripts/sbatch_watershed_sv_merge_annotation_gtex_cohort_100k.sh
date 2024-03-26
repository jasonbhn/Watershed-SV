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
#SBATCH --partition=defq
#SBATCH --mem=100G
#SBATCH --signal=USR2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
module load anaconda
conda activate
conda activate rare_variants
module load bcftools
module load samtools
module load bedtools

python3.9 scripts/executable_scripts/combine_all_annotations.py \
--vcf input/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes new-protein-lincRNA-100k-7.0/intermediates/pipeline_input_genotypes.tsv \
--genes new-protein-lincRNA-100k-7.0/intermediates/genes.bed \
--gene-sv new-protein-lincRNA-100k-7.0/intermediates/gene_sv.100000.bed \
--annotation-dir new-protein-lincRNA-100k-7.0/intermediates/ \
--outfile new-protein-lincRNA-100k-7.0/combined_annotation_pre_merge_gene_level_noimpute.csv \
--expressions proteinCoding_lincRNA/intermediate/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.z3.txt \
--expression-field MedZ \
--expression-id-field Ind \
--maf-mode upload \
--maf-file new-protein-lincRNA-100k-7.0/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file new-protein-lincRNA-100k-7.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 100000

python3.9 scripts/executable_scripts/combine_all_annotations.py \
--vcf input/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes new-protein-lincRNA-100k-7.0/intermediates/pipeline_input_genotypes.tsv \
--genes new-protein-lincRNA-100k-7.0/intermediates/genes.bed \
--gene-sv new-protein-lincRNA-100k-7.0/intermediates/gene_sv.100000.bed \
--annotation-dir new-protein-lincRNA-100k-7.0/intermediates/ \
--outfile new-protein-lincRNA-100k-7.0/combined_annotation_pre_merge_noimpute.csv \
--expressions proteinCoding_lincRNA/intermediate/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.z3.txt \
--expression-field MedZ \
--expression-id-field Ind \
--maf-mode upload \
--maf-file new-protein-lincRNA-100k-7.0/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file new-protein-lincRNA-100k-7.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000