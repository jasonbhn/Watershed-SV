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
#SBATCH --time=3:0:0
#SBATCH --mem=200G
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
mamba activate rare_variants
module load bcftools
module load samtools
module load bedtools
bash scripts/executable_scripts/generate_annotations_ABC.sh \
--pipeline smallset \
--input-vcf CMG_input/SV_calls_Parliament2_SVTools_concat_w_INS/CMG_Parliament2_SVTools_jasmine_wINS.cohort.sorted.MEIasDEL.noIllegal.reIDed.final.rare.vcf \
--filters "PASS" \
--flank 100000 \
--medz-cutoff 3 \
--rareness 0.01 \
--outdir MuscularDystrophy-protein-lincRNA-100k-all-position-aware-3.0/ \
--genome-bound-file input/gtex_human_reference.genome \
--transcript-origin input/gencode.v26.GRCh38.genes.gtf \
--gencode-genes input/gencode.v26.GRCh38.genes.gtf \
--genome-build hg38 \
--vep-cache-dir vep_cache \
--filter-ethnicity False 

