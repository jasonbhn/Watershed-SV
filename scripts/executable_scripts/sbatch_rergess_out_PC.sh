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

python3.9 scripts/executable_scripts/regress_out_PC.py \
--regression-data /scratch16/abattle4/bohan/Watershed-SV/input/CMG_MD_PCA_study_regression_data_with_age_disease_status_120PCs_no_bad_samples.csv \
--num-PC-remove $1 \
--out-parquet-dir /data/abattle4/bohan/Muscular_dystrophy_expression/CMG_MD_PCA_study_regress_top_$1_PCs_no_bad_samples/
