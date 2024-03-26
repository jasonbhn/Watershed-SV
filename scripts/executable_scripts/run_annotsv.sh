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
#SBATCH --time=5:0:0
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
# input is a csv of input, output, HPO. 
module load anaconda
conda activate
conda activate rare_variants

while read -r input output hpo_terms
do
$ANNOTSV/bin/AnnotSV -SVinputFile $input \
-outputFile $output \
-hpo $hpo_terms \
-svtBEDcol 4 \
-samplesidBEDcol 5
done < "$1"
