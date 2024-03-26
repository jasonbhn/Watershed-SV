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
cp input/gtex.lumpy.gs.melt.high_conf.vcf.gz /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz
cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/input/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.z3.standard.txt /tmp/TMM_corrected_zscores_tissue.csv
cp -r new-protein-lincRNA-100k-rare-position-aware-2.0 /tmp/

python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/ \
--outfile new-protein-lincRNA-100k-rare-position-aware-2.0/combined_annotation_pre_merge_gene_level_noimpute.medZ.csv \
--expressions /tmp/TMM_corrected_zscores_tissue.csv \
--expression-field MedZ \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 100000 \
--zscore-threshold 3 \
--minimum-support-tissue-count 1

python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/ \
--outfile new-protein-lincRNA-100k-rare-position-aware-2.0/combined_annotation_pre_merge_noimpute.medZ.csv \
--expressions /tmp/TMM_corrected_zscores_tissue.csv \
--expression-field MedZ \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-position-aware-2.0/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000 \
--zscore-threshold 3 \
--minimum-support-tissue-count 1
