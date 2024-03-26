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
#SBATCH --time=5:0:0
#SBATCH --mem=128G
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
#module load ensembl-vep/109.3
mamba activate rare_variants
module load bcftools
module load samtools
module load bedtools
bash scripts/executable_scripts/generate_annotations.sh \
--pipeline smallset \
--input-vcf CMG_input/SV_calls_Parliament2_SVTools_concat_w_INS/CMG_Parliament2_SVTools_jasmine_wINS.cohort.sorted.MEIasDEL.noIllegal.reIDed.final.rare.vcf \
--filters "PASS" \
--flank 10000 \
--medz-cutoff 3 \
--rareness 0.01 \
--outdir MuscularDystrophy-protein-lincRNA-10k-final-version/ \
--genome-bound-file input/gtex_human_reference.genome \
--transcript-origin input/gencode.v26.GRCh38.genes.gtf \
--gencode-genes input/gencode.v26.GRCh38.genes.gtf \
--genome-build hg38 \
--filter-ethnicity False \
--vep-cache-dir /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/vep_cache
