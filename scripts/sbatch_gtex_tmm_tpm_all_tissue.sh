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
#SBATCH --mem=64G
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
mamba activate r_env
#declare -a tissues=('Muscle_Skeletal' 'Whole_Blood' 'Skin_Sun_Exposed_Lower_leg' 
#'Adipose_Subcutaneous' 'Artery_Tibial' 'Thyroid' 'Nerve_Tibial' 
#'Skin_Not_Sun_Exposed_Suprapubic' 'Lung' 'Esophagus_Mucosa' 
#'Adipose_Visceral_Omentum' 'Esophagus_Muscularis' 'Cells_Cultured_fibroblasts' 
#'Breast_Mammary_Tissue' 'Heart_Left_Ventricle' 'Artery_Aorta' 
#'Heart_Atrial_Appendage' 'Colon_Transverse' 'Esophagus_Gastroesophageal_Junction' 
#'Colon_Sigmoid' 'Testis' 'Stomach' 'Pancreas' 'Pituitary' 'Adrenal_Gland' 
#'Brain_Cortex' 'Brain_Caudate_basal_ganglia' 'Brain_Nucleus_accumbens_basal_ganglia' 
#'Prostate' 'Brain_Cerebellum' 'Spleen' 'Artery_Coronary' 'Liver' 'Brain_Cerebellar_Hemisphere' 
#'Brain_Frontal_Cortex_BA9' 'Brain_Putamen_basal_ganglia' 'Brain_Hypothalamus' 'Brain_Hippocampus' 
#'Small_Intestine_Terminal_Ileum' 'Ovary' 'Brain_Anterior_cingulate_cortex_BA24' 
#'Cells_EBV-transformed_lymphocytes' 'Minor_Salivary_Gland' 'Brain_Spinal_cord_cervical_c-1' 'Vagina' 
#'Brain_Amygdala' 'Uterus' 'Brain_Substantia_nigra' 'Kidney_Cortex')
declare -a tissues=('Muscle_Skeletal' 'Whole_Blood' 'Cells_Cultured_fibroblasts')
for tissue in "${tissues[@]}"
do
 Rscript scripts/executable_scripts/Expression_outlier_calling_in_one_cohort.R -t "$tissue" -o input/TMM_corrected_tissue_GTEx_no_batch/
 echo "$tissue"
done
