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
#SBATCH --time=0:15:0
#SBATCH --mem=300G
#SBATCH --signal=USR2
#SBATCH --nodes=1
#SBATCH --account=smontgom
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
module load anaconda
eval "$(conda shell.bash hook)"
conda activate rare_variants

python3.9 scripts/executable_scripts/enrichment_by_rare_SV.py --input-file-list UDN_expression_list.txt --outfile UDN-protein-lincRNA-10k-10.0/outlier_enrichment_for_rareSV.parquet --sv-annot-dir UDN-protein-lincRNA-10k-10.0/ --maf-file UDN-protein-lincRNA-10k-10.0/intermediates/custom_UDN_maf.tsv --flank 10000 --no-use-gtex --af-colname af --sample-colname UDNID --expression-colname zscore
