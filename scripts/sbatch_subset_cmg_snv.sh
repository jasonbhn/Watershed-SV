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
#SBATCH --time=3:0:0
#SBATCH --partition=bigmem
#SBATCH -A mschatz1_bigmem
#SBATCH --mem=200G
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
module load bcftools
module load samtools
module load bedtools

bcftools index -t /scratch16/abattle4/lab_data/muscular_dystrophy/GATK_SNV_calls_postVSQR/CMG.filterPASS.new.vcf.gz
#bcftools view -i 'FILTER="PASS"' -s 143BN_BB,167BZ_SP,179CI_GG,251DW_SD,253DY_HA,255EA_GC,65T_CR,79Z_CP,B09-25.1,B10-02.1,B11-25.1,B11-48.1,B12-21.1,B13-07.1,B13-15,B13-52.1,B14-07,B14-117.1,B14-130.1,B14-48.1,BON_B09-27.1,BON_B12-74.1,BON_B14-71.2,BON_UC219.1,UC223.1,UC305.1,UC316.1,UC84.1 -Oz -o /scratch16/abattle4/lab_data/muscular_dystrophy/GATK_SNV_calls_postVSQR/CMG.filterPASS.new.vcf.gz /scratch16/abattle4/lab_data/muscular_dystrophy/GATK_SNV_calls_postVSQR/CMG_padded.vcf.gz
