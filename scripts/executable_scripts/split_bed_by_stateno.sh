#!/usrbin/env bash

set -o nounset -o pipefail -o errexit

indir=$1
outdir=$2
if [ ! -d "$2" ] 
then 
    echo "$outdir not here, let's create it."
    mkdir $outdir
fi
for tissue_epigenome in $indir/*.bed
do    
    for state in {1..25}
    do  
        #echo "processing $tissue_epigenome state number $state"
        awk -v s=$state '$4==s {print $0}' $tissue_epigenome > $outdir/$(basename $tissue_epigenome .bed).$state.bed
        if [ -s "$outdir/$(basename $tissue_epigenome .bed).$state.bed" ] 
        then
            echo "---"
        else
            echo "$outdir/$(basename $tissue_epigenome .bed).$state.bed is empty"
        fi
    done
done
        
