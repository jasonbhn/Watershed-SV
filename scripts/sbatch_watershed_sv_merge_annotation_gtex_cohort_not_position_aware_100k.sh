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
cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/input/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.z3.standard.txt /tmp/
cp -r new-protein-lincRNA-100k-rare-baseline /tmp/
cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/input/TMM_corrected_zscores_Muscle_Skeletal.tsv /tmp/
cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/input/TMM_corrected_zscores_Whole_Blood.tsv /tmp/
cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/input/TMM_corrected_zscores_Cells_Cultured_fibroblasts.tsv /tmp/

 

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_gene_level_noimpute.medZ.csv \
--expressions /tmp/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.z3.standard.txt \
--expression-field MedZ \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 100000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_noimpute.medZ.csv \
--expressions /tmp/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.z3.standard.txt \
--expression-field MedZ \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000

# muscle
python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_gene_level_noimpute.Muscle_Skeletal.csv \
--expressions /tmp/TMM_corrected_zscores_Muscle_Skeletal.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 100000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_noimpute.Muscle_Skeletal.csv \
--expressions /tmp/TMM_corrected_zscores_Muscle_Skeletal.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_gene_level_noimpute.Whole_Blood.csv \
--expressions /tmp/TMM_corrected_zscores_Whole_Blood.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 100000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_noimpute.Whole_Blood.csv \
--expressions /tmp/TMM_corrected_zscores_Whole_Blood.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_gene_level_noimpute.Cells_Cultured_fibroblasts.csv \
--expressions /tmp/TMM_corrected_zscores_Cells_Cultured_fibroblasts.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene \
--remove-control-genes \
--flank 100000

python3.10 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /tmp/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--genotypes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/genes.bed \
--gene-sv /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/gene_sv.100000.bed \
--annotation-dir /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/ \
--outfile new-protein-lincRNA-100k-rare-baseline/combined_annotation_pre_merge_noimpute.Cells_Cultured_fibroblasts.csv \
--expressions /tmp/TMM_corrected_zscores_Cells_Cultured_fibroblasts.tsv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode upload \
--maf-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_gtex_maf.tsv \
--length-mode extract \
--length-field SVLEN \
--CN-mode upload \
--CN-file /tmp/new-protein-lincRNA-100k-rare-baseline/intermediates/custom_CN_frame.tsv \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000
