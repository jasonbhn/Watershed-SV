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
#SBATCH --time=0:30:0
#SBATCH --account=smontgom
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
module load anaconda
eval "$(conda shell.bash hook)"
conda activate rare_variants

python3.9 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN365359/933309-UDN365359-P.SV_reheadered.svafotate.ofp95.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/ \
--outfile UDN-protein-lincRNA-10k-UDN365359-1.0/combined_annotation_pre_merge_gene_Blood_custom_maf.csv \
--expressions UDN_dataset/Blood.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.reIDed.txt \
--expression-field zscore \
--expression-id-field UDNID \
--maf-mode upload \
--maf-file UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/custom_MAF.tsv \
--length-mode upload \
--length-file UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/custom_length.tsv \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000

python3.9 scripts/executable_scripts/combine_all_annotations_polar.py \
--vcf /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN365359/933309-UDN365359-P.SV_reheadered.svafotate.ofp95.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/ \
--outfile UDN-protein-lincRNA-10k-UDN365359-1.0/combined_annotation_pre_merge_gene_Fibroblast_custom_maf.csv \
--expressions UDN_dataset/Fibroblast.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.reIDed.txt \
--expression-field zscore \
--expression-id-field UDNID \
--maf-mode upload \
--maf-file UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/custom_MAF.tsv \
--length-mode upload \
--length-file UDN-protein-lincRNA-10k-UDN365359-1.0/intermediates/custom_length.tsv \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000
