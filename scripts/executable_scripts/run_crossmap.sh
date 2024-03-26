#!/bin/bash
#SBATCH --job-name=sd_0.4_0.1
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --mem=100G
ml gcc/5.5.0
ml python/3.7-anaconda-2019.03
source ~/.bashrc
conda activate
conda activate rare_variants
CrossMap.py bigwig /work-zfs/abattle4/bohan/SV-watershed/Annotations/hg19ToHg38.over.chain.gz /work-zfs/abattle4/bohan/SV-watershed/starting_annotations/LINSIGHT.bw /work-zfs/abattle4/bohan/SV-watershed/starting_annotations/LINSIGHT38.test
conda deactivate
