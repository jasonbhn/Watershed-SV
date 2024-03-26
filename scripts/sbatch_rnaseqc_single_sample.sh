#!/bin/bash
# ---------------------------------------------------
#  INPUT ENVIRONMENT VARIABLES
# ---------------------------------------------------
#SBATCH --job-name=WtsdSV
#SBATCH --time=16:0:0
#SBATCH --partition=batch
#SBATCH -A smontgom
#SBATCH --mem=32G
#SBATCH --signal=USR2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------

set -e
# step 1. parse arguments: use getopt
SHORT=g:,s:,b:,o:,h
LONG=gene-gtf:,sample-id:,bam:,outdir:,help
OPTS=$(getopt -a -n sv_annotations --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

while :
do
  case "$1" in
    -g | --gene-gtf)
      gene_gtf="$2"
      shift 2
      ;;
    -s | --sample-id )
      sample_id="$2"
      shift 2
      ;;
    -b | --bam )
      bam="$2"
      shift 2
      ;;
    -o | --outdir )
      outdir="$2"
      shift 2
      ;;
    -h | --help)
      "RNASeQC quantify genes"
      exit 2
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      ;;
  esac
done
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD433_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD433
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD434_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD434
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD435_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD435
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD436_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD436
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD437_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD437
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD438_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD438
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD439_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD439
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD440_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD440
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD441_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD441
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD442_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD442
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD443_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD443
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD444_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD444
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD445_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD445
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD446_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD446
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD447_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD447
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD448_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD448
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD452_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD452
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD454_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD454
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD462_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD462
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD463_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD463
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD464_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD464
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD465_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD465
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD466_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD466
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD467_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD467
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD469_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD469
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD470_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD470
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD471_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD471
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD533_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD533
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD534_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD534
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD535_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD535
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD536_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD536
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD537_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD537
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD538_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD538
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD539_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD539
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD540_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD540
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD541_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD541
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD542_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD542
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD543_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD543
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD544_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD544
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD545_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD545
rnaseqc.v2.4.2.linux input/gencode.v35.primary_assembly.gene.gtf /oak/stanford/groups/smontgom/shared/UDN/Old/PreprocessingHG38Primary/BAMs/RD546_star_hg38_Aligned.sortedByCoord.sorted_opticalDup12000_dedupOptical_minMQ255.bam UDN_dataset/UDN_RNASeQC/ -s RD546
