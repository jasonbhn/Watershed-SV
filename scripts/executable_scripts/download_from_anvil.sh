#!/bin/bash
set -e
# step 1. parse arguments: use getopt
SHORT=l:,o:,h
LONG=urilist:,outdir:,help
OPTS=$(getopt -a -n sv_annotations --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

while :
do
  case "$1" in
    -l | --urilist )
      urilist="$2"
      shift 2
      ;;
    -o | --outdir )
      outdir="$2"
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

# download file from anvil
while read p; do
  echo "$p"
  gsutil -m cp $p $outdir
done <$urilist