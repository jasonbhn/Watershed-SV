#!/bin/bash
# ---------------------------------------------------
#  INPUT ENVIRONMENT VARIABLES
# ---------------------------------------------------
#SBATCH --job-name=WtsdSV
#SBATCH --time=1:0:0
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
ml rsem/1.3.3
set -e
# step 1. parse arguments: use getopt
SHORT=r:,s:,b:,o:,h
LONG=rsem-reference:,sample-id:,bam:,outdir:,help
OPTS=$(getopt -a -n sv_annotations --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

while :
do
  case "$1" in
    -r | --rsem-reference)
      rsem_ref="$2"
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
      "RSEM quantify genes"
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

echo "starting rsem quantification..."
rsem-calculate-expression -p 4 \
--paired-end \
--keep-intermediate-files \
--estimate-rspd \
--strandedness none \
--seed 12345 \
--bam $bam \
$rsem_ref \
$outdir
echo "done!"
