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
#SBATCH --time=1:0:0
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --account=smontgom
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
module load anaconda
eval "$(conda shell.bash hook)"
conda activate rare_variants
module load bcftools
module load samtools
module load bedtools

python3.9 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf CMG_input/SV_calls_Parliament2_processed/$1.survivor.qual.filterPASS.2callers.rare.vcf \
--genotypes MuscularDystrophy-protein-lincRNA-10k-single-sample-1.0/$1/intermediates/pipeline_input_genotypes.tsv \
--genes MuscularDystrophy-protein-lincRNA-10k-single-sample-1.0/$1/intermediates/genes.bed \
--gene-sv MuscularDystrophy-protein-lincRNA-10k-single-sample-1.0/$1/intermediates/gene_sv.10000.bed \
--annotation-dir MuscularDystrophy-protein-lincRNA-10k-single-sample-1.0/$1/intermediates/ \
--outfile MuscularDystrophy-protein-lincRNA-10k-single-sample-1.0/$1/combined_annotation_pre_merge_gene_sv_new_expression.csv \
--expressions CMG_input/CMG_83_PCs_only.expression.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode extract \
--maf-field Max_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000

#python3.9 scripts/executable_scripts/train_test_predict_split_annotation.py \
#--training new-protein-lincRNA-100k-7.0/combined_annotation_pre_merge_gene_level_noimpute.csv \
#--testing-list MuscularDystrophy-protein-lincRNA-10k-single-sample/CMG_annotations_list.csv \
#--mode predict \
#--min-af-impute-mode infer \
#--out-prefix MuscularDystrophy-protein-lincRNA-10k-single-sample/CMG_parliament_lumpy

#sbatch test_drive_predict.sh ../MuscularDystrophy-protein-lincRNA-10k-single-sample/CMG_parliament_lumpy_training_data.tsv ../MuscularDystrophy-protein-lincRNA-10k-single-sample/CMG_parliament_lumpy_testing_data.tsv ../MuscularDystrophy-protein-lincRNA-10k-single-sample/CMG_parliament_lumpy_prediction 10000 0.0027 2.56 1 RIVER
