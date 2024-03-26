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
SHORT=r:,p:,o:,h
LONG=result-dir:,prefix:,outdir:,help
OPTS=$(getopt -a -n sv_annotations --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

while :
do
  case "$1" in
    -r | --result-dir)
      result_dir="$2"
      shift 2
      ;;
    -p | --prefix )
      prefix="$2"
      shift 2
      ;;
    -o | --outdir )
      outdir="$2"
      shift 2
      ;;
    -h | --help)
      "RNASeQC aggregate"
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

module load anaconda
eval "$(conda shell.bash hook)"
conda activate rare_variants

python3.9 -m rnaseqc aggregate --parquet -o $outdir $result_dir $prefix
