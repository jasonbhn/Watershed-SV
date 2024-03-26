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
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
eval "$(conda shell.bash hook)"
conda activate rare_variants

python3.10 scripts/executable_scripts/train_test_predict_split_annotation.py --training new-protein-lincRNA-10k-9.0/combined_annotation_pre_merge_gene_level_noimpute.csv --mode predict --out-prefix UDN-protein-lincRNA-10k-ADRC-new/Fibro_annotations_custom_maf --min-af-impute-mode infer --testings UDN-protein-lincRNA-10k-ADRC-new/combined_annotation_pre_merge_gene_Fibroblast_custom_maf.csv
python3.10 scripts/executable_scripts/train_test_predict_split_annotation.py --training new-protein-lincRNA-10k-9.0/combined_annotation_pre_merge_gene_level_noimpute.csv --mode predict --out-prefix UDN-protein-lincRNA-10k-ADRC-new/Blood_annotations_custom_maf --min-af-impute-mode infer --testings UDN-protein-lincRNA-10k-ADRC-new/combined_annotation_pre_merge_gene_Blood_custom_maf.csv
