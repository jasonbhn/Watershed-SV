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
--vcf /scratch16/abattle4/bohan/Watershed-SV/rare_disease_validation/SV_calls/cmg_muscle_disease_27_indiv.all.filterPASS.svafotate.vcf.gz \
--genotypes MuscularDystrophy-protein-lincRNA-100k-4.0/intermediates/pipeline_input_genotypes.tsv \
--genes MuscularDystrophy-protein-lincRNA-100k-4.0/intermediates/genes.bed \
--gene-sv MuscularDystrophy-protein-lincRNA-100k-4.0/intermediates/gene_sv.100000.bed \
--annotation-dir MuscularDystrophy-protein-lincRNA-100k-4.0/intermediates/ \
--outfile MuscularDystrophy-protein-lincRNA-100k-4.0/combined_annotation_pre_merge_gene.csv \
--expressions /scratch16/abattle4/lab_data/muscular_dystrophy/cmg_expression.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode extract \
--maf-field Max_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file MuscularDystrophy-protein-lincRNA-100k-4.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000
