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
#SBATCH --mem=180G
#SBATCH --signal=USR2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
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
--vcf /scratch16/abattle4/bohan/Watershed-SV/rare_disease_validation/SV_calls/UDN/UDN.cohort_combined.GRCh38.consensus_SVs.jasmine_irisRefined.major_chr_seq_resolved.paragraph_1kgp_maf.reheader.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-8.0/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-8.0/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-8.0/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-8.0/intermediates/ \
--outfile UDN-protein-lincRNA-10k-8.0/combined_annotation_pre_merge_gene_Blood.csv \
--expressions UDN_dataset/Blood.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.reIDed.txt \
--expression-field zscore \
--expression-id-field UDNID \
--maf-mode extract \
--maf-field 1KGP_PARAGRAPH_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000

python3.9 scripts/executable_scripts/combine_all_annotations.py \
--vcf /scratch16/abattle4/bohan/Watershed-SV/rare_disease_validation/SV_calls/UDN/UDN.cohort_combined.GRCh38.consensus_SVs.jasmine_irisRefined.major_chr_seq_resolved.paragraph_1kgp_maf.reheader.vcf.gz \
--genotypes UDN-protein-lincRNA-10k-8.0/intermediates/pipeline_input_genotypes.tsv \
--genes UDN-protein-lincRNA-10k-8.0/intermediates/genes.bed \
--gene-sv UDN-protein-lincRNA-10k-8.0/intermediates/gene_sv.10000.bed \
--annotation-dir UDN-protein-lincRNA-10k-8.0/intermediates/ \
--outfile UDN-protein-lincRNA-10k-8.0/combined_annotation_pre_merge_gene_Fibroblast.csv \
--expressions UDN_dataset/Fibroblast.hg38.sorted_opticalDupdedupOptical_minMQ255_rsem.genes.zscores.reIDed.txt \
--expression-field zscore \
--expression-id-field UDNID \
--maf-mode extract \
--maf-field 1KGP_PARAGRAPH_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--flank 10000
