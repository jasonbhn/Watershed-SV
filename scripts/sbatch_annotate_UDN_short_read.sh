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
#SBATCH --time=4:0:0
#SBATCH --account=smontgom
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
module load anaconda
eval "$(conda shell.bash hook)"
conda activate rare_variants
svafotate annotate -v ../UDN.COHORT_MERGED.hg38.short_read_wgs.SVs.AF_annotated.recompressed.reheader.correctheader.vcf.gz -o ../UDN.COHORT_MERGED.hg38.short_read_wgs.SVs.AF_annotated.recompressed.reheader.correctheader.svafotate090.vcf.gz -b ../SVAFotate_core_SV_popAFs.GRCh38.bed.gz -f 0.90 -a best
svafotate annotate -v ../UDN.COHORT_MERGED.hg38.short_read_wgs.SVs.AF_annotated.recompressed.reheader.correctheader.vcf.gz -o ../UDN.COHORT_MERGED.hg38.short_read_wgs.SVs.AF_annotated.recompressed.reheader.correctheader.svafotate095.vcf.gz -b ../SVAFotate_core_SV_popAFs.GRCh38.bed.gz -f 0.95 -a best
