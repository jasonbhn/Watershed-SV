# Watershed-SV
Watershed-SV extends [Watershed](https://github.com/BennyStrobes/Watershed) to model the impact of rare SVs (DUP, DEL, DUP-CNV, DEL-CNV, INV, INS) on nearby gene expressions outlier. For running Watershed model, please refer to Watershed GitHub. This repository contains:
1. the pipeline and associated scripts used for generating structural variations (rare and common) annotations with respect to nearby genes.
2. the scripts for generating expression outliers. 
3. the approach to merge annotations AND the expression outliers, and finally format data into desired format for [evaluate_watershed.R](https://github.com/BennyStrobes/Watershed/blob/master/evaluate_watershed.R) and [predict_watershed.R](https://github.com/BennyStrobes/Watershed/blob/master/predict_watershed.R).

# Evaluating Watershed-SV model against WGS-only model

# Using Watershed-SV model learned on a reference dataset to prioritize gene-SV pairs in an independent test dataset. 

