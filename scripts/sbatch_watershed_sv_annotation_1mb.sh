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
#SBATCH --mem=300G
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
--pipeline population \
--input-vcf input/gtex.lumpy.gs.melt.high_conf.vcf.gz \
--filters "MATCH_1KGP . PASS" \
--flank 2000000 \
--medz-cutoff 3 \
--rareness 0.01 \
--liftover-bed input/gtex.lumpy.gs.melt.high_conf.hg38.bed \
--outdir new-protein-lincRNA-1000k-rare-1.0/ \
--genome-bound-file input/gtex_human_reference.genome \
--transcript-origin input/gencode.v26.GRCh38.genes.gtf \
--gencode-genes input/gencode.v26.GRCh38.genes.gtf \
--genome-build hg38 \
--vep-cache-dir /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/vep_cache \
--metadata sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt \
--filter-ethnicity True \
--filter-rare True
