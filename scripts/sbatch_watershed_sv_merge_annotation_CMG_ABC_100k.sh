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
#SBATCH --time=2:0:0
#SBATCH --account=smontgom
#SBATCH --mem=128G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base

# Set bash options for verbose output and to fail immediately on errors or if variables are undefined.
mamba activate rare_variants


cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/CMG_input/SV_calls_Parliament2_SVTools_concat_w_INS/CMG_Parliament2_SVTools_jasmine_wINS.cohort.sorted.MEIasDEL.noIllegal.reIDed.final.rare.vcf /tmp/ 
cp UDN_dataset/finalized_UDN_expression/combined_Muscle_top_60PC_RIN_SEX_TMM_TPM.nodup.noGlobalOutliers.zscore.tsv /tmp/
cp input/CMG_hybrid_maf.wRepeat.tsv /tmp/
cp -r MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/ /tmp/ 

python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/CMG_Parliament2_SVTools_jasmine_wINS.cohort.sorted.MEIasDEL.noIllegal.reIDed.final.rare.vcf \
--genotypes /tmp/MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/intermediates/genes.bed \
--gene-sv /tmp/MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/ \
--outfile MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/combined_annotation_pre_merge_custom_maf.no_control_gene.tmm_corrected.z3.csv \
--expressions /tmp/combined_Muscle_top_60PC_RIN_SEX_TMM_TPM.nodup.noGlobalOutliers.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/CMG_hybrid_maf.wRepeat.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000 \
--zscore-threshold 3 \
--minimum-support-tissue-count 1


python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/CMG_Parliament2_SVTools_jasmine_wINS.cohort.sorted.MEIasDEL.noIllegal.reIDed.final.rare.vcf \
--genotypes /tmp/MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/intermediates/genes.bed \
--gene-sv /tmp/MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/ \
--outfile MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/combined_annotation_pre_merge_custom_maf_gene_level.no_control_gene.tmm_corrected.z3.csv \
--expressions /tmp/combined_Muscle_top_60PC_RIN_SEX_TMM_TPM.nodup.noGlobalOutliers.zscore.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/CMG_hybrid_maf.wRepeat.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene \
--remove-control-genes \
--flank 100000 \
--zscore-threshold 3 \
--minimum-support-tissue-count 1

