#!/bin/bash
set -e
# step 1. parse arguments: use getopt
SHORT=v:,o:,s:,h
LONG=vcf:,outdir:,subjectID:,help
OPTS=$(getopt -a -n sv_annotations --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

while :
do
  case "$1" in
    -v | --vcf )
      vcf="$2"
      shift 2
      ;;
    -o | --outdir )
      outdir="$2"
      shift 2
      ;;
    -s | --subjectID )
      subjectID="$2"
      shift 2
      ;;
    -h | --help)
      echo "Extract genomic annotations for SVs v1"
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
mkdir -p $outdir
vcfdir="$(dirname "${vcf}")"
vcfbase="$(basename "${vcf}" .vcf)"
bcftools view -Oz -i 'FILTER="PASS"&&ABS(SVLEN)<10000000' $vcf > $vcfdir/$vcfbase.filterPASS.preCallerFilter.vcf.gz
python3.9 /scratch16/abattle4/bohan/Watershed-SV/scripts/executable_scripts/filterByCaller.py --vcf $vcfdir/$vcfbase.filterPASS.preCallerFilter.vcf.gz --outfile $vcfdir/$vcfbase.filterPASS.2callers.vcf.gz --subjectID $subjectID
svafotate annotate --minf 0.9 -a best -v $vcfdir/$vcfbase.filterPASS.2callers.vcf.gz -o $vcfdir/$vcfbase.filterPASS.2callers.svafotate.vcf -b /scratch16/abattle4/bohan/Watershed-SV/rare_disease_validation/SV_calls/SVAFotate_core_SV_popAFs.GRCh38.bed.gz
#get bed file containing variant loc, variant genotype
bcftools query -i \
'FILTER="PASS"&Max_AF!=0' -f \
'%CHROM\t%POS0\t%END\t%ID\t%INFO/SubjectID\t[%GT\t%SP]\t%SVTYPE\t%SVLEN\t%Max_AF\n' $vcfdir/$vcfbase.filterPASS.2callers.svafotate.vcf > \
$outdir/$vcfbase.filterPASS.2callers.svafotate.known.bed
bcftools query -i \
'FILTER="PASS"&Max_AF==0' -f \
'%CHROM\t%POS0\t%END\t%ID\t%INFO/SubjectID\t[%GT\t%SP]\t%SVTYPE\t%SVLEN\t%Max_AF\n' $vcfdir/$vcfbase.filterPASS.2callers.svafotate.vcf > \
$outdir/$vcfbase.filterPASS.2callers.svafotate.putative_novel.bed
bcftools view -Ov -i 'FILTER="PASS"&Max_AF==0' -o $outdir/$vcfbase.filterPASS.2callers.svafotate.putative_novel.vcf $vcfdir/$vcfbase.filterPASS.2callers.svafotate.vcf
python3.9 /scratch16/abattle4/bohan/Watershed-SV/scripts/executable_scripts/convert_vcf_to_paragraph_input.py \
--vcf $outdir/$vcfbase.filterPASS.2callers.svafotate.putative_novel.vcf --outfile $outdir/$vcfbase.filterPASS.2callers.svafotate.putative_novel.paragraph.vcf



