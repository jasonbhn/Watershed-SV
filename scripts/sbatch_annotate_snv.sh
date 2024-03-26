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
#SBATCH --time=6:0:0
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --account=smontgom
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
# Run scripts to enable activating conda environments
. "${HOME}/micromamba/etc/profile.d/conda.sh"
. "${HOME}/micromamba/etc/profile.d/mamba.sh"
mamba activate base
module load ensembl-vep/109.3
mamba activate rare_variants
module load bcftools
module load samtools
module load bedtools


bash scripts/executable_scripts/extract_snv_vep_annotations.sh /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN022194_Clair3_vcf/UDN022194_Fibro.GRCh38.clair3.SNP_indel.whatshap_PHASED.vcf.gz 
