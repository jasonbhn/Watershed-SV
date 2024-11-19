#!/usr/bin/env bash

set -o nounset -o pipefail -o errexit

inputDir=$1
outputDir=$2
if [ ! -d "$outputDir" ] 
then
    echo "Directory $outputDir DOES NOT exists. Creating it" 
    mkdir $outputDir
fi
for raw_tad in $inputDir/*.txt; 
do 
    awk '{OFS="\t";if ($3 > $2 && $3 >= 0 && $2 >= 0 ) print $1,$2,$2+1,NR-1; else if ($3 < $2 && $3 >= 0 && $2 >= 0) print $1,$3,$3+1,NR-1;}' $raw_tad | sort -k1,1 -k2,2n >  $outputDir/$(basename "$raw_tad" .txt).start.bed
    awk '{OFS="\t";if ($3 > $2 && $3 >= 0 && $2 >= 0 ) print $1,$3,$3+1,NR-1; else if ($3 < $2 && $3 >= 0 && $2 >= 0) print $1,$2,$2+1,NR-1;}' $raw_tad | sort -k1,1 -k2,2n >  $outputDir/$(basename "$raw_tad" .txt).end.bed
done
