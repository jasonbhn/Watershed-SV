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
which python
bash scripts/executable_scripts/generate_annotations_ABC.sh \
--pipeline smallset \
--input-vcf ../UDN_MergedSV_vcfs/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.noBNDTRA.vcf.gz \
--filters "PASS" \
--flank 100000 \
--medz-cutoff 3 \
--rareness 0.01 \
--outdir UDN-protein-lincRNA-100k-ADRC-all-position-aware-2.0/ \
--genome-bound-file input/gtex_human_reference.genome \
--transcript-origin input/gencode.v35.primary_assembly.annotation.gtf \
--gencode-genes input/gencode.v35.primary_assembly.annotation.gtf \
--genome-build hg38 \
--vep-cache-dir vep_cache \
--filter-ethnicity False
