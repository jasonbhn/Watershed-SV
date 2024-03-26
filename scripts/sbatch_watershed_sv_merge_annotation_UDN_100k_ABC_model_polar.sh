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
#SBATCH --time=2:0:0
#SBATCH --account=smontgom
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
eval "$(conda shell.bash hook)"
conda activate rare_variants

cp /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN_MergedSV_vcfs/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.noBNDTRA.vcf.gz /tmp/ 
cp input/UDN_LRS_ABC_udn_expression.tmm_corrected.csv /tmp/
cp -r UDN-protein-lincRNA-100k-ADRC-all-position-aware-1.0/ /tmp/ 
echo "here"
python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
--vcf /tmp/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.noBNDTRA.vcf.gz \
--genotypes /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-1.0/intermediates/pipeline_input_genotypes.tsv \
--genes /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-1.0/intermediates/genes.bed \
--gene-sv /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-1.0/intermediates/gene_sv_slop.100000.bed \
--annotation-dir /tmp/UDN-protein-lincRNA-100k-ADRC-all-position-aware-1.0/ \
--outfile UDN-protein-lincRNA-100k-ADRC-all-position-aware-1.0/combined_annotation_pre_merge_custom_maf.no_control_gene.tmm_corrected.no_vep.csv \
--expressions /tmp/UDN_LRS_ABC_udn_expression.tmm_corrected.csv \
--expression-field normalized_resid \
--expression-id-field Ind \
--maf-mode extract \
--maf-field ALNG_AF \
--length-mode extract \
--length-field SVLEN \
--CN-mode extract \
--collapse-mode gene-sv \
--remove-control-genes \
--flank 100000

python3.10 scripts/executable_scripts/train_test_predict_split_annotation.py --training new-protein-lincRNA-100k-rare-position-aware-1.0/combined_annotation_pre_merge_gene_level_noimpute.tissue.tmm_corrected.z3.csv --testings UDN-protein-lincRNA-100k-ADRC-all-position-aware-1.0/combined_annotation_pre_merge_custom_maf.no_control_gene.tmm_corrected.no_vep.csv --mode predict --min-af-impute-mode infer --out-prefix UDN-protein-lincRNA-100k-ADRC-all-position-aware-1.0/tissue_specific_annotations_no_vep_sv_pred
