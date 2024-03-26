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
#SBATCH --partition=batch
#SBATCH --account=smontgom
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
eval "$(conda shell.bash hook)"
conda activate rare_variants

python3.9 scripts/executable_scripts/regress_out_PC_no_age.py \
--regression-data /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/Blood.expression.hg38.gtex.udn.noDup.tsv \
--num-PC-remove 60 \
--out-parquet-dir /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/UDN_blood_regress_top_60_PCs_Sex/

python3.9 scripts/executable_scripts/regress_out_PC_no_age.py \
--regression-data /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/Fibroblast.expression.hg38.gtex.udn.noDup.tsv \
--num-PC-remove 60 \
--out-parquet-dir /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/UDN_fibroblast_regress_top_60_PCs_Sex/
