# Watershed-SV [Construction in progress]
Watershed-SV extends [Watershed](https://github.com/BennyStrobes/Watershed) to model the impact of rare SVs (DUP, DEL, DUP-CNV, DEL-CNV, INV, INS) on nearby gene expressions outlier. For running Watershed model, please refer to Watershed GitHub. This repository contains:
1. the pipeline and associated scripts used for generating structural variations (rare and common) annotations with respect to nearby genes.
2. the scripts for generating expression outliers. 
3. the approach to merge annotations AND the expression outliers, and finally format data into desired format for [evaluate_watershed.R](https://github.com/BennyStrobes/Watershed/blob/master/evaluate_watershed.R) and [predict_watershed.R](https://github.com/BennyStrobes/Watershed/blob/master/predict_watershed.R).

## Instructions for collecting annotations. 
Watershed-SV pipeline currently uses bash script for simplicity. The key scripts are `scripts/executable_scripts/generate_annotations.sh` and `scripts/executable_scripts/generate_annotations_ABC.sh` depending on whether you want to run region-agnostic model (10kb model in the paper) or region-aware model (100kb model in the paper). 
To replicate the environment for collecting annotations, see: `WatershedSV.yml`. 
### Parameter breakdown: 
1. `-p | --pipeline`: which pipeline to use, select from `population`, `smallset`.
If data is of sufficient size, ie > 100, select population, allowing for option --filter-ethnicity, --filter rare. 
Otherwise, select smallset. 
2. `-v | --input-vcf`: input vcf file, it has to have at least 1 sample column. We only consider SVTYPEs: DUP, DEL, DUP_CNV, DEL_CNV, CNV, INS, INV. 
3. `-f | --filters`: if variant record in vcf have these filters, keep for further analysis. 
4. `-k | --flank`: how much flanking up and downstream of genes to consider. usually use 100000, 10000.  
5. `-r | --rareness`: rareness, if --filter-rare == True, then this is the MAF threshold to set to filter for rare variants, 0.01 recommended or lower.  
6. `-l | --liftover-bed`: if you have a crossmap/liftover SV coordinate you want to use, ie, if VCF is in older build, you lifted over coordinates to HG38, then provide the bed file in addition to original VCF to convert coordinates. 
7. `-o | --outdir`: output directory name for annotations
8. `-b | --genome-bound-file`: a file depicting the chromosome/contig name, start and end coordinates.
9. `-g | --gencode-genes`: gencode transcript model file.
10. `-c | --vep-cache-dir`: vep_cache_dir for running vep annotations. we recommend setting up vep offline to run our pipeline smoothly.
11. `-a | --metadata`: metadata file for filtering ethnicity. In our case, training data is GTEx, we used GTEx metadata file from dbGaP.
12. `-e | --filter-ethnicity`: filter by ethnicity? GTEx relic, True to only train on EUR individuals.
13. `-i | --filter-rare`: filter rare variants if using `population` model 

## Evaluating Watershed-SV model against WGS-only model

## Using Watershed-SV model learned on a reference dataset to prioritize gene-SV pairs in an independent test dataset. 

