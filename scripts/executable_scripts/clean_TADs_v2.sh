#!/bin/bash

inputDir=$1
outputDir=$2
if [ ! -d "$outputDir" ] 
then
    echo "Directory $outputDir DOES NOT exists. Creating it" 
    mkdir $outputDir
fi
for raw_tad in $inputDir/*.txt; 
do 
    awk \
    '{
        OFS="\t";
        if($3 > $2)
            print $1,$2,$3,NR-1;
        else
            print $1,$3,$2,NR-1;}' \
    $raw_tad | sort -k1,1 -k2,2n >  $outputDir/$(basename "$raw_tad" .txt).bed
    #awk '$3 > $2{OFS="\t";print $1,$3-1,$3,NR-1}' $raw_tad | sort -k1,1 -k2,2n >  $outputDir/$(basename "$raw_tad" .txt).end.bed
done
