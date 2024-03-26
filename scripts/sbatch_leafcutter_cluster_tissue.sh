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
#SBATCH --time=12:0:0
#SBATCH --partition=defq
#SBATCH --account=mschatz1
#SBATCH --mem=32G
#SBATCH --signal=USR2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bni1@jhu.edu
#SBATCH --output=Watershed-SV.job.%j.out
#SBATCH --erro=Watershed-SV.job.%j.err
# ---------------------------------------------------
module load anaconda
conda activate
conda activate rare_variants
set -e
# step 1. parse arguments: use getopt
SHORT=j:,w:,h
LONG=junclist:,workdir:,help
OPTS=$(getopt -a -n sv_annotations --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

while :
do
  case "$1" in
    -j | --junclist )
      junclist="$2"
      shift 2
      ;;
    -w | --workdir )
      workdir="$2"
      shift 2
      ;;
    -h | --help)
      "Extract genomic annotations for SVs v1"
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
# initial clustering. 
python /scratch16/abattle4/bohan/Watershed-SV/scripts/executable_scripts/leafcutter_cluster_regtools.py \
-j $junclist \
-m 50 \
-r $workdir \
-o tmp \
-l 500000
# now remove all of the bad exon exon junctions and recluster. 
python3.9 /scratch16/abattle4/bohan/Watershed-SV/scripts/executable_scripts/filter_splice_clusters.py \
--cluster $workdir/tmp_perind_numers.counts.gz \
--filtered-junclist $workdir/tmp_junclist.txt \
--junclist $junclist
# post junction filter clustering. 
python /scratch16/abattle4/bohan/Watershed-SV/scripts/executable_scripts/leafcutter_cluster_regtools.py \
-j $workdir/tmp_junclist.txt \
-m 50 \
-r $workdir \
-o tmp_filtered \
-l 500000
