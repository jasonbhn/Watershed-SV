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
#SBATCH --time=10:0:0
#SBATCH --mem=32G
#SBATCH --account=smontgom
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
echo $1
bash scripts/executable_scripts/generate_annotations.sh \
--pipeline smallset \
--input-vcf CMG_input/SV_calls_Parliament2_processed/$1.survivor.qual.filterPASS.2callers.rare.vcf \
--filters "PASS" \
--flank 10000 \
--medz-cutoff 3 \
--rareness 0.01 \
--liftover-bed input/gtex.lumpy.gs.melt.high_conf.hg38.bed \
--outdir MuscularDystrophy-protein-lincRNA-10k-single-sample-1.0/$1 \
--genome-bound-file input/gtex_human_reference.genome \
--transcript-origin input/gencode.v26.GRCh38.genes.gtf \
--gencode-genes input/gencode.v26.GRCh38.genes.gtf \
--genome-build hg38 \
--vep-cache-dir vep_cache \
--filter-ethnicity False
