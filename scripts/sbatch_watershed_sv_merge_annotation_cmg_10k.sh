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
#SBATCH --mem=100G
#SBATCH --account=smontgom
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
module load anaconda
eval "$(conda shell.bash hook)"
conda activate rare_variants


python3.9 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf CMG_input/CMG_Parliament2_SVTools_jasmine_wINS.cohort.sorted.MEIasDEL.RepeatMasked.noIllegal.final.reIDed.vcf \
--genotypes MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/pipeline_input_genotypes.tsv \
--genes MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/genes.bed \
--gene-sv MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/gene_sv.10000.bed \
--annotation-dir MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/ \
--outfile MuscularDystrophy-protein-lincRNA-10k-6.0/combined_annotation_pre_merge_gene_sv_expnDiseaseSexPC_strictAF.csv \
--expressions CMG_input/CMG_81_PCs_disease_sex.expression.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/custom_CMG_MAF_localAC_stringent.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 10000
python3.9 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf CMG_input/cmg_muscle_disease_27_indiv.all.filterPASS.svafotate.vcf.gz \
--genotypes MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/pipeline_input_genotypes.tsv \
--genes MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/genes.bed \
--gene-sv MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/gene_sv.10000.bed \
--annotation-dir MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/ \
--outfile MuscularDystrophy-protein-lincRNA-10k-6.0/combined_annotation_pre_merge_gene_sv_expnDiseaseSexPC_strictAF_gene_level.csv \
--expressions CMG_input/CMG_81_PCs_disease_sex.expression.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/custom_CMG_MAF_localAC_stringent.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file MuscularDystrophy-protein-lincRNA-10k-6.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 10000
