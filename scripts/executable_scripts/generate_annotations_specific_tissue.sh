#!/bin/bash
set -e
# step 1. parse arguments: use getopt
SHORT=v:,k:,m:,o:,b:,t:,g:,s:,d:,c:,h
LONG=gene-sv:,flank:,medz-cutoff:,outdir:\
,genome-bound-file:,transcript-origin:,gencode-genes:,gene-list:,genome-build:,vep-cache-dir:,help
OPTS=$(getopt -a -n sv_annotations --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

while :
do
  case "$1" in
    -v | --gene-sv )
      flank="$2"
      shift 2
      ;;
    -k | --flank )
      flank="$2"
      shift 2
      ;;
    -m | --medz-cutoff )
      medz_cutoff="$2"
      shift 2
      ;;
    -o | --outdir )
      outdir="$2"
      shift 2
      ;;
    -b | --genome-bound-file )
      genome_bound_file="$2"
      shift 2
      ;;
    -t | --transcript-origin )
      transcript_origin="$2"
      shift 2
      ;;
    -g | --gencode-genes )
      gencode_genes="$2"
      shift 2
      ;;
    -s | --gene-list )
      gene_list="$2"
      shift 2
      ;;
    -d | --genome-build )
      genome_build="$2"
      shift 2
      ;;
    -c | --vep-cache-dir )
      vep_cache_dir="$2"
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
# make directories for outputs.
mkdir -p $outdir
mkdir -p $outdir/intermediates

# gene annotations processing
python scripts/executable_scripts/extract_gene_exec.py \
--gencode-annotations $gencode_genes \
--gene-list $gene_list \
--out-gene-bed ${outdir}/intermediates/genes.bed \
--out-exon-bed ${outdir}/intermediates/exons.bed \
--out-gene-tss ${outdir}/intermediates/gene_tss.tsv \
--out-gene-tes ${outdir}/intermediates/gene_tes.tsv
echo 'Done extracting gene annotations'
# step 2. get non personal SV annotation, annotations based on SV coordinates and potentially content. 
gene_sv=/scratch16/abattle4/bohan/Watershed-SV/new-protein-lincRNA/intermediates/gene_sv_cmg.10000.bed

# sv_to_exon:
bedtools intersect -wa -wb -a ${outdir}/intermediates/exons.bed -b ${gene_sv} | 
awk 'BEGIN{{print "SV\tGene\texon"}};$4==$11{{OFS="\t";print $9,$11,$5}}' > ${outdir}/intermediates/exon_sv.${flank}.tsv
echo 'Done extracting exon sv overlaps'

# sv_to_gene_remap:
remap_crm=input/remap2020_crm_macs2_hg38_v1_0.bed.gz
zcat ${remap_crm} | awk '{OFS="\t";print $1,$2,$3,$5}' | 
bedtools intersect -wa -wb -a stdin -b ${gene_sv} | 
awk 'BEGIN{print "SV\tGene\tremap_crm_score"};{OFS="\t";print $8,$10,$4}' > ${outdir}/intermediates/remap_crm_sv.dist.${flank}.tsv
echo 'Done extracting gene sv remap overlap'
# TODO add conditional block for tissue specific annotations

# sv_to_gene_dist:
python scripts/executable_scripts/sv_to_gene_dist.py \
--gene-sv ${gene_sv} \
--in-gene-tss ${outdir}/intermediates/gene_tss.tsv \
--in-gene-tes ${outdir}/intermediates/gene_tes.tsv \
--out-gene-sv-dist ${outdir}/intermediates/SV_dist_to_gene.dist.${flank}.tsv
echo 'Done extracting gene sv dist'

# sv_to_gene_cpg_shell:
cpg=input/cpgIslandExt.txt
# out format: 'chrom','start','end','SV','Gene','chrom','start','end','cpg'
awk '{{FS="\t";OFS="\t";print $2,$3,$4,$10}}' ${cpg} > ${outdir}/intermediates/cpgtmp.bed
bedtools intersect -wa -wb -a ${gene_sv} -b ${outdir}/intermediates/cpgtmp.bed > ${outdir}/intermediates/cpg_by_genes_SV.dist.${flank}.bed

# sv_to_gene_cpg_py:
python scripts/executable_scripts/sv_to_gene_cpg.py --gene-sv-cpg ${outdir}/intermediates/cpg_by_genes_SV.dist.${flank}.bed --out-gene-sv-cpg ${outdir}/intermediates/sv_to_gene_cpg.dist.${flank}.tsv
echo 'Done extracting gene sv cpg'

# sv_to_gene_bw_scores_py:
# gc content
gc=input/gc5Base.bw
python scripts/executable_scripts/sv_to_gene_bw_scores.py \
--gene-sv ${gene_sv} \
--in-bigwig ${gc} \
--bigwig-name "mean_GC_content" \
--stat-method "mean" \
--score-upper-limit 100 \
--score-lower-limit 0 \
--out-gene-sv-score ${outdir}/intermediates/GC_by_genes_SV.dist.${flank}.tsv
echo 'Done extracting gene sv gc'

# linsight scores
linsight_file=input/LINSIGHT_hg38.bw
python scripts/executable_scripts/sv_to_gene_bw_scores.py \
--gene-sv ${gene_sv} \
--in-bigwig ${linsight_file} \
--bigwig-name "mean_LINSIGHT" \
--stat-method "mean" \
--score-upper-limit 1 \
--score-lower-limit 0 \
--out-gene-sv-score ${outdir}/intermediates/LINSIGHT_by_genes_SV.dist.${flank}.tsv
echo 'Done extracting gene sv linsight'

# phastCON 20 scores
phcon20_file=input/hg38.phastCons20way.bw
python scripts/executable_scripts/sv_to_gene_bw_scores.py \
--gene-sv ${gene_sv} \
--in-bigwig ${phcon20_file} \
--bigwig-name "mean_phastCON" \
--stat-method "mean" \
--score-upper-limit 1 \
--score-lower-limit 0 \
--out-gene-sv-score ${outdir}/intermediates/PhastCON20_by_genes_SV.dist.${flank}.tsv
echo 'Done extracting gene sv phastcon'

# CADD scores
cadd_file=input/CADD_GRCh38-v1.6.bw # use mean top 10
python scripts/executable_scripts/sv_to_gene_bw_scores.py \
--gene-sv ${gene_sv} \
--in-bigwig ${phcon20_file} \
--bigwig-name "top10mean_CADD" \
--stat-method "top10_mean" \
--score-upper-limit 1 \
--score-lower-limit 0 \
--out-gene-sv-score ${outdir}/intermediates/CADD_by_genes_SV.dist.${flank}.tsv
echo 'Done extracting gene sv CADD'


# sv_to_gene_TADs:
TADs_dir=input/TAD # hg38 from starting annotations
bash scripts/executable_scripts/clean_TADs_v2.sh ${TADs_dir} ${outdir}/intermediates/TAD_cleaned
bedtools slop -i ${gene_sv} -g ${genome_bound_file} -b 5000 | 
bedtools intersect -wa -wb -a stdin -b ${outdir}/intermediates/TAD_cleaned/* -filenames | 
sort -k4,4 -k6,6 | 
bedtools groupby -i stdin -g 4,6 -c 7 -o count_distinct | 
awk 'BEGIN{{print "SV\tGene\tnum_TADs"}};{{OFS="\t";print}}' > ${outdir}/intermediates/TAD_boundary_by_genes_SV.dist.${flank}.tsv
echo 'Done extracting gene sv TAD'

# merge_enhancers:
enhancers=input/matrix_hs.csv
primary_cells=input/hs_primary.txt

python scripts/executable_scripts/merge_enhancers.py \
--enhancers ${enhancers} \
--primary-cell-list ${primary_cells} \
--out-merged-enhancers ${outdir}/intermediates/primary_cells_collpased_enhancers.bed
# sv_to_gene_enhancers:
bedtools intersect -wa -wb -a ${outdir}/intermediates/primary_cells_collpased_enhancers.bed -b ${gene_sv} | 
awk 'BEGIN{{print "SV\tGene\tnum_enhancers_cell_types"}};{{OFS="\t";print $8,$10,$4}}' > ${outdir}/intermediates/enhancers_by_genes_SV.dist.${flank}.tsv
echo 'Done extracting gene sv enhancer'

# process_roadmaps:
roadmap_dir=input/roadmap_all_gtex
scripts/executable_scripts/split_bed_by_stateno.sh ${roadmap_dir} ${outdir}/intermediates/processed_roadmaps
# sv_to_gene_roadmaps
for i in {1..25}
do
scripts/executable_scripts/sv_to_gene_roadmap.sh ${gene_sv} ${outdir}/intermediates/processed_roadmaps ${outdir}/intermediates/roadmap_multitissue_sv_to_gene.${i}.tsv ${i}
done

# combine_sv_to_gene_roadmaps:
python scripts/executable_scripts/combine_roadmaps.py \
--gene-sv-roadmap-dir ${outdir}/intermediates \
--out-combined-roadmap ${outdir}/intermediates/combined_roadmaps.dist.${flank}.tsv
echo 'Done extracting gene sv roadmaps'

# vep annotations
python scripts/executable_scripts/prep_vep_input.py ${gene_sv} ${outdir}/intermediates/vep_input.${flank}.bed
scripts/executable_scripts/extract_sv_vep_annotations.sh /home/bni1/vep.sif \
${outdir}/intermediates/vep_input.${flank}.bed \
${vep_cache_dir} \
${outdir}/intermediates/tmp.tsv \
${outdir}/intermediates/sv_to_gene_vep.${flank}.tsv \
${flank}


# now ready to combine


