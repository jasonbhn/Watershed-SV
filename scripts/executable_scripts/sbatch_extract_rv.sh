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
#SBATCH --time=3:0:0
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
python scripts/executable_scripts/extract_rare_variants.py \
--vcf /scratch16/abattle4/bohan/Watershed-SV/input/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--lifted-coord /scratch16/abattle4/bohan/Watershed-SV/input/gtex.lumpy.gs.melt.high_conf.hg38.bed \
--extract-genotype \
--infer-rareness \
--genotype-filters PASS . MATCH_1KGP \
--out-annotsv /scratch16/abattle4/bohan/Watershed-SV/new-protein-lincRNA/intermediates/vep_input.tsv \
--out-generic /scratch16/abattle4/bohan/Watershed-SV/new-protein-lincRNA/intermediates/pipeline_input.bed \
--out-maf /scratch16/abattle4/bohan/Watershed-SV/new-protein-lincRNA/intermediates/pipeline_maf.tsv \
--out-genotype /scratch16/abattle4/bohan/Watershed-SV/new-protein-lincRNA/intermediates/pipeline_input_genotypes.tsv


