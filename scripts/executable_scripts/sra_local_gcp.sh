#!/bin/bash

Help()
{
   # Display Help
   echo "Download fastq from SRA, compress them to fastq.gz, then upload to Terra."
   echo
   echo "Syntax: sra_local_gcp.sh  [-s|m|o|b|h]"
   echo "options:"
   echo "-s     --sra-accession: accession list file. one SRRXXXXX per line"
   echo "-m     --multiprocess-limit: max number of parallel at the same time"
   echo "-o     --outdir: output directory"
   echo "-n     --ngc: path to ngc key"
   echo "-b     --bucket: which bucket directory to upload to"
   echo "-h     --help"
   echo
}
# parse args
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--sra-accession) sra_list_file="$2"; shift ;;
        -m|--multiprocess-limit) m="$2"; shift ;;
        -o|--outdir) outdir="$2"; shift ;;
        -n|--ngc) ngc="$2"; shift ;;
        -b|--bucket) bucket="$2"; shift ;;
        -h|--help) Help; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done
# make outdir if not exist
mkdir -p $outdir
echo "starting downlink from sra to local"
# start download using pipe to parallel
while read p; do
if [[ -f $outdir/"$p"_1.fastq ]] 
then
    echo "$p exist!"
else
    fasterq-dump --threads 4 --ngc $ngc -O $outdir $p
fi
done < $sra_list_file
#cat $sra_list_file | parallel -j $m fasterq-dump --threads 4 --ngc $ngc -O $outdir
echo "done downloading, now compress"
gzip $outdir/*.fastq
echo "done downloading,compressing, now start uplink with gcp"
gsutil -m cp $outdir $bucket
echo "transfer complete"

