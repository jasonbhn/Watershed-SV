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
#SBATCH --time=1:0:0
#SBATCH --account=smontgom
#SBATCH --mem=100G
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
module load bcftools
module load samtools
module load bedtools
module load ensembl-vep/109.3

bash scripts/executable_scripts/generate_annotations.sh \
--pipeline smallset \
--input-vcf ../UDN365359/933309-UDN365359-P.SV_reheadered.vcf.gz \
--filters "PASS" \
--flank 10000 \
--medz-cutoff 3 \
--rareness 0.01 \
--outdir UDN-protein-lincRNA-10k-UDN365359-1.0/ \
--genome-bound-file input/gtex_human_reference.genome \
--transcript-origin input/gencode.v42.basic.annotation.gtf \
--gencode-genes input/gencode.v42.basic.annotation.gtf \
--genome-build hg38 \
--vep-cache-dir vep_cache \
--filter-ethnicity False
