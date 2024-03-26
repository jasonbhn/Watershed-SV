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
#SBATCH --time=12:0:0
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

bash /scratch16/abattle4/bohan/Watershed-SV/scripts/executable_scripts/generate_annotations.sh \
--pipeline smallset \
--input-vcf /scratch16/abattle4/bohan/Watershed-SV/rare_disease_validation/SV_calls/UDN/UDN.cohort_combined.GRCh38.consensus_SVs.jasmine_irisRefined.major_chr_seq_resolved.paragraph_1kgp_maf.reheader.vcf.gz \
--filters "PASS" \
--flank 100000 \
--medz-cutoff 3 \
--rareness 0.01 \
--outdir /scratch16/abattle4/bohan/Watershed-SV/UDN-protein-lincRNA-100k-8.0/ \
--genome-bound-file /scratch16/abattle4/bohan/Watershed-SV/input/gtex_human_reference.genome \
--transcript-origin /scratch16/abattle4/bohan/Watershed-SV/input/gencode.v42.basic.annotation.gtf \
--gencode-genes /scratch16/abattle4/bohan/Watershed-SV/input/gencode.v42.basic.annotation.gtf \
--genome-build hg38 \
--vep-cache-dir /scratch16/abattle4/bohan/vep_data \
--filter-ethnicity False
