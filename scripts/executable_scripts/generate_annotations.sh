#!/usr/bin/env bash

set -o nounset -o pipefail -o errexit
# step 1. parse arguments: use getopt
SHORT=p:,v:,f:,k:,m:,r:,l:,o:,b:,t:,g:,d:,c:,a:,e:,i:,h
LONG=pipeline:,input-vcf:,filters:,flank:,medz-cutoff:,rareness:,liftover-bed:,outdir:\
,genome-bound-file:,transcript-origin:,gencode-genes:,genome-build:,vep-cache-dir:,metadata:,filter-ethnicity:,filter-rare:,help
OPTS=$(getopt -a -n sv_annotations --options $SHORT --longoptions $LONG -- "$@")

eval set -- "$OPTS"

while :
do
  case "$1" in
    -p | --pipeline )
      pipeline="$2"
      shift 2
      ;;
    -v | --input-vcf )
      input_vcf="$2"
      shift 2
      ;;
    -f | --filters )
      filters="$2"
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
    -r | --rareness )
      rareness="$2"
      shift 2
      ;;
    -l | --liftover-bed )
      liftover_bed="$2"
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
    -d | --genome-build )
      genome_build="$2"
      shift 2
      ;;
    -c | --vep-cache-dir )
      vep_cache_dir="$2"
      shift 2
      ;;
    -a | --metadata )
      metadata="$2"
      shift 2
      ;;
    -e | --filter-ethnicity )
      filter_ethnicity="$2"
      shift 2
      ;;
    -i | --filter-rare )
      filter_rare="$2"
      shift 2
      ;;
    -h | --help )
      echo "Extract genomic annotations for SVs v1"
      exit 0
      ;;
    -- )
      shift;
      break
      ;;
    * )
      echo "Unexpected option: $1"
      exit 1
      ;;
  esac
done
# make directories for outputs.
mkdir -p $outdir
mkdir -p $outdir/intermediates
# step 1. extract variants bed files for annotations based on pipeline mode
# available modes: population, smallset, sample
if [ "$pipeline" == "population" ]; then
  if [ ! -f "${outdir}/intermediates/pipeline_input.bed" ]; then
    if [ "$filter_ethnicity" == "True" ]; then
      if [ "$filter_rare" == "True" ]; then 
        extract_rare_variants \
        --vcf $input_vcf \
        --lifted-coord $liftover_bed \
        --extract-genotype \
        --infer-rareness \
        --filter-ethnicity \
        --metadata $metadata \
        --genotype-filters $filters \
        --out-annotsv ${outdir}/intermediates/vep_input.tsv \
        --out-generic ${outdir}/intermediates/pipeline_input.bed \
        --out-maf ${outdir}/intermediates/pipeline_maf.tsv \
        --out-genotype ${outdir}/intermediates/pipeline_input_genotypes.tsv
      else 
        extract_rare_variants \
        --vcf $input_vcf \
        --lifted-coord $liftover_bed \
        --extract-genotype \
        --filter-ethnicity \
        --metadata $metadata \
        --genotype-filters $filters \
        --out-annotsv ${outdir}/intermediates/vep_input.tsv \
        --out-generic ${outdir}/intermediates/pipeline_input.bed \
        --out-maf ${outdir}/intermediates/pipeline_maf.tsv \
        --out-genotype ${outdir}/intermediates/pipeline_input_genotypes.tsv
      fi
    else
      if [ "$filter_rare" == "True" ]; then 
        extract_rare_variants \
        --vcf $input_vcf \
        --lifted-coord $liftover_bed \
        --extract-genotype \
        --infer-rareness \
        --genotype-filters $filters \
        --out-annotsv ${outdir}/intermediates/vep_input.tsv \
        --out-generic ${outdir}/intermediates/pipeline_input.bed \
        --out-maf ${outdir}/intermediates/pipeline_maf.tsv \
        --out-genotype ${outdir}/intermediates/pipeline_input_genotypes.tsv
      else
        extract_rare_variants \
        --vcf $input_vcf \
        --lifted-coord $liftover_bed \
        --extract-genotype \
        --genotype-filters $filters \
        --out-annotsv ${outdir}/intermediates/vep_input.tsv \
        --out-generic ${outdir}/intermediates/pipeline_input.bed \
        --out-maf ${outdir}/intermediates/pipeline_maf.tsv \
        --out-genotype ${outdir}/intermediates/pipeline_input_genotypes.tsv
      fi
    fi
  fi
elif [ "$pipeline" == "smallset" ]; then
  if [ ! -f "${outdir}/intermediates/pipeline_input.bed" ]; then
  extract_rare_variants \
  --vcf $input_vcf \
  --extract-genotype \
  --genotype-filters ${filters} \
  --out-annotsv ${outdir}/intermediates/vep_input.tsv \
  --out-generic ${outdir}/intermediates/pipeline_input.bed \
  --out-genotype ${outdir}/intermediates/pipeline_input_genotypes.tsv
  fi
  
fi


# gene annotations processing
if [ ! -f "${outdir}/intermediates/genes.bed" ]; then
extract_gene_exec \
--gencode-annotations $gencode_genes \
--out-gene-bed ${outdir}/intermediates/genes.bed \
--out-exon-bed ${outdir}/intermediates/exons.bed \
--out-gene-tss ${outdir}/intermediates/gene_tss.tsv \
--out-gene-tes ${outdir}/intermediates/gene_tes.tsv \
--genome-bound ${genome_bound_file}
fi
echo 'Done extracting gene annotations'
# step 2. get non personal SV annotation, annotations based on SV coordinates and potentially content. 

# sv_to_gene processing
if [ ! -f "${outdir}/intermediates/gene_sv.${flank}.bed" ]; then
bedtools slop -g ${genome_bound_file} -i ${outdir}/intermediates/genes.bed -b ${flank} | 
bedtools intersect -a ${outdir}/intermediates/pipeline_input.bed -b stdin -wb | 
awk '{{OFS="\t";print $1,$2,$3,$4,$5,$9}}' > ${outdir}/intermediates/gene_sv.${flank}.bed
fi
echo 'Done extracting gene sv overlaps'


# sv_to_exon:
if [ ! -f "${outdir}/intermediates/exon_sv.${flank}.tsv" ]; then
bedtools intersect -wo -a ${outdir}/intermediates/exons.bed -b ${outdir}/intermediates/gene_sv.${flank}.bed | 
awk '$4==$14{OFS="\t";print $4,$12,$5,$6,$7,$8,$15}' | 
sort -k1,1 -k2,2 -k3,3n > ${outdir}/intermediates/exon_sv.${flank}.unprocessed_info.tsv
extract_SV_exon_info \
--input ${outdir}/intermediates/exon_sv.${flank}.unprocessed_info.tsv \
--output ${outdir}/intermediates/exon_sv.${flank}.tsv
fi
echo 'Done extracting exon sv overlaps'

# sv_to_gene_remap:
if [ ! -f "${outdir}/intermediates/remap_crm_sv.dist.${flank}.tsv" ]; then
remap_crm=input/remap2020_crm_macs2_hg38_v1_0.bed.gz
zcat ${remap_crm} | awk '{OFS="\t";print $1,$2,$3,$5}' | 
bedtools intersect -wa -wb -a stdin -b ${outdir}/intermediates/gene_sv.${flank}.bed | 
awk '{OFS="\t";print $8,$10,$4}' | 
sort -k1,1 -k2,2 | 
bedtools groupby -i stdin -g 1,2 -c 3 -o max | 
awk 'BEGIN{print "SV\tGene\tremap_crm_score"};{OFS="\t";print $0}' \
> ${outdir}/intermediates/remap_crm_sv.dist.${flank}.tsv
fi
echo 'Done extracting gene sv remap overlap'
# TODO add conditional block for tissue specific annotations

# sv_to_gene_dist:
if [ ! -f "${outdir}/intermediates/SV_dist_to_gene.dist.${flank}.tsv" ]; then
sv_to_gene_dist \
--flank ${flank} \
--gene ${outdir}/intermediates/genes.bed \
--gene-sv ${outdir}/intermediates/gene_sv.${flank}.bed \
--in-gene-tss ${outdir}/intermediates/gene_tss.tsv \
--in-gene-tes ${outdir}/intermediates/gene_tes.tsv \
--out-gene-sv-dist ${outdir}/intermediates/SV_dist_to_gene.dist.${flank}.tsv
fi
echo 'Done extracting gene sv dist'

# ABC model 
if [ ! -f "${outdir}/intermediates/SV_dist_to_gene.dist.${flank}.tsv" ]; then
ABC_enhancers=input/ABC_enhancers.merged.All.sorted.11col.region.hg38.bed
bedtools intersect -wa -wb -a ${outdir}/intermediates/gene_sv.${flank}.bed -b ${ABC_enhancers} |
awk 'BEGIN{print "SV\tGene\tis_ABC_SV"}; split($6,a,".") {OFS="\t";if (a[1]==$17) print $4,a[1],1}' > ${outdir}/intermediates/sv_to_gene_ABC.${flank}.tsv
fi
echo 'Done extracting gene sv abc annotations'
# sv_to_gene_cpg_shell:
cpg=input/cpgIslandExt.txt
if [ ! -f "${outdir}/intermediates/sv_to_gene_cpg.dist.${flank}.tsv" ]; then
# out format: 'chrom','start','end','SV','Gene','chrom','start','end','cpg'
awk '{{FS="\t";OFS="\t";print $2,$3,$4,$10}}' ${cpg} > ${outdir}/intermediates/cpgtmp.bed
bedtools intersect -wa -wb -a ${outdir}/intermediates/gene_sv.${flank}.bed -b ${outdir}/intermediates/cpgtmp.bed > ${outdir}/intermediates/cpg_by_genes_SV.dist.${flank}.bed

# sv_to_gene_cpg_py:
sv_to_gene_cpg --gene-sv-cpg ${outdir}/intermediates/cpg_by_genes_SV.dist.${flank}.bed --out-gene-sv-cpg ${outdir}/intermediates/sv_to_gene_cpg.dist.${flank}.tsv
fi
echo 'Done extracting gene sv cpg'

# sv_to_gene_bw_scores_py:
# gc content
if [ ! -f "${outdir}/intermediates/GC_by_genes_SV.dist.${flank}.tsv" ]; then
gc=input/gc5Base.bw
sv_to_gene_bw_scores \
--gene-sv ${outdir}/intermediates/gene_sv.${flank}.bed \
--in-bigwig ${gc} \
--bigwig-name "mean_GC_content" \
--stat-method "mean" \
--score-upper-limit 100 \
--score-lower-limit 0 \
--out-gene-sv-score ${outdir}/intermediates/GC_by_genes_SV.dist.${flank}.tsv
fi
echo 'Done extracting gene sv gc'

# linsight scores
if [ ! -f "${outdir}/intermediates/LINSIGHT_by_genes_SV.dist.${flank}.tsv" ]; then
linsight_file=input/LINSIGHT_hg38.bw
sv_to_gene_bw_scores \
--gene-sv ${outdir}/intermediates/gene_sv.${flank}.bed \
--in-bigwig ${linsight_file} \
--bigwig-name "top10_LINSIGHT" \
--stat-method "top10_mean" \
--score-upper-limit 1 \
--score-lower-limit 0 \
--out-gene-sv-score ${outdir}/intermediates/LINSIGHT_by_genes_SV.dist.${flank}.tsv
fi
echo 'Done extracting gene sv linsight'

# phastCON 20 scores
if [ ! -f "${outdir}/intermediates/PhastCON20_by_genes_SV.dist.${flank}.tsv" ]; then
phcon20_file=input/hg38.phastCons20way.bw
sv_to_gene_bw_scores \
--gene-sv ${outdir}/intermediates/gene_sv.${flank}.bed \
--in-bigwig ${phcon20_file} \
--bigwig-name "top10_phastCON" \
--stat-method "top10_mean" \
--score-upper-limit 1 \
--score-lower-limit 0 \
--out-gene-sv-score ${outdir}/intermediates/PhastCON20_by_genes_SV.dist.${flank}.tsv
fi
echo 'Done extracting gene sv phastcon'

# CADD scores
if [ ! -f "${outdir}/intermediates/CADD_by_genes_SV.dist.${flank}.tsv" ]; then
cadd_file=input/CADD_GRCh38-v1.6.bw # use mean top 10
sv_to_gene_bw_scores \
--gene-sv ${outdir}/intermediates/gene_sv.${flank}.bed \
--in-bigwig ${cadd_file} \
--bigwig-name "top10_CADD" \
--stat-method "top10_mean" \
--score-upper-limit 99 \
--score-lower-limit 0 \
--out-gene-sv-score ${outdir}/intermediates/CADD_by_genes_SV.dist.${flank}.tsv
fi
echo 'Done extracting gene sv CADD'


# sv_to_gene_TADs:
if [ ! -f "${outdir}/intermediates/TAD_boundary_by_genes_SV.dist.${flank}.tsv" ]; then
TADs_dir=input/TAD # hg38 from starting annotations
clean_TADs ${TADs_dir} ${outdir}/intermediates/TAD_cleaned
bedtools slop -i ${outdir}/intermediates/gene_sv.${flank}.bed -g ${genome_bound_file} -b 5000 > ${outdir}/intermediates/TAD_5000_flank_gene_SV.bed
bedtools intersect -wa -wb -a ${outdir}/intermediates/TAD_5000_flank_gene_SV.bed -b ${outdir}/intermediates/TAD_cleaned/*  -filenames | 
sort -k4,4 -k6,6 | 
bedtools groupby -i stdin -g 4,6 -c 7 -o count_distinct | 
awk 'BEGIN{{print "SV\tGene\tnum_TADs"}};{{OFS="\t";print}}' > ${outdir}/intermediates/TAD_boundary_by_genes_SV.dist.${flank}.tsv
echo 'Done extracting gene sv TAD'
fi
# merge_enhancers:
enhancers=input/matrix_hs.csv
primary_cells=input/hs_primary.txt
if [ ! -f "${outdir}/intermediates/enhancers_by_genes_SV.dist.${flank}.tsv" ]; then
merge_enhancers \
--enhancers ${enhancers} \
--primary-cell-list ${primary_cells} \
--out-merged-enhancers ${outdir}/intermediates/primary_cells_collpased_enhancers.bed
# sv_to_gene_enhancers:
bedtools intersect -wa -wb -a ${outdir}/intermediates/primary_cells_collpased_enhancers.bed -b ${outdir}/intermediates/gene_sv.${flank}.bed | 
awk '{OFS="\t";print $8,$10,$4}' | 
sort -k1,1 -k2,2 | 
bedtools groupby -i stdin -g 1,2 -c 3 -o max | 
awk 'BEGIN{print "SV\tGene\tnum_enhancers_cell_types"};{OFS="\t";print $0}' \
> ${outdir}/intermediates/enhancers_by_genes_SV.dist.${flank}.tsv
fi
echo 'Done extracting gene sv enhancer'

# process_roadmaps:
if [ ! -f "${outdir}/intermediates/combined_roadmaps.dist.${flank}.tsv" ]; then
roadmap_dir=input/roadmap_all_gtex
split_bed_by_stateno ${roadmap_dir} ${outdir}/intermediates/processed_roadmaps
# sv_to_gene_roadmaps
for i in {1..25}
do
sv_to_gene_roadmap ${outdir}/intermediates/gene_sv.${flank}.bed ${outdir}/intermediates/processed_roadmaps ${outdir}/intermediates/roadmap_multitissue_sv_to_gene.${i}.tsv ${i}
done

# combine_sv_to_gene_roadmaps:
combine_roadmaps \
--gene-sv-roadmap-dir ${outdir}/intermediates \
--out-combined-roadmap ${outdir}/intermediates/combined_roadmaps.dist.${flank}.tsv
fi
echo 'Done extracting gene sv roadmaps'

# vep annotations
if [ ! -f "${outdir}/intermediates/sv_to_gene_vep.${flank}.tsv" ]; then
prep_vep_input \
${outdir}/intermediates/gene_sv.${flank}.bed \
${outdir}/intermediates/vep_input.${flank}.bed

run_extract_sv_vep_annotations \
${outdir}/intermediates/vep_input.${flank}.bed \
${vep_cache_dir} \
${outdir}/intermediates/tmp.tsv \
${outdir}/intermediates/sv_to_gene_vep.${flank}.tsv \
fi
echo 'Done extracting gene sv vep annotations'
