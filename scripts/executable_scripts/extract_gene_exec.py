#!/usr/bin/env python3

import pandas as pd
import polars as pl
import argparse
import pyranges as pr
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gene related annotation processing')
    parser.add_argument('--gencode-annotations',type=str,required=True,metavar='[gencode annotation file]',
                    help='gencode annotations to get gene information in gtf')
    parser.add_argument('--out-gene-bed',type=str,required=True,metavar='[input for SV-gene pair annotations]',
                    help='path to input file to be created for SV-gene annotations')
    parser.add_argument('--out-exon-bed',type=str,required=True,metavar='[input for SV-exon pair annotations]',
                    help='path to input file to be created for SV-exon annotations')
    parser.add_argument('--out-gene-tss',type=str,required=True,metavar='[input for SV-tss pair annotations]',
                    help='path to input file to be created for SV-tss annotations')
    parser.add_argument('--out-gene-tes',type=str,required=True,metavar='[input for SV-tes pair annotations]',
                    help='path to input file to be created for SV-tes annotations')
    parser.add_argument('--genome-bound',type=str,required=True,metavar='[input for chrom bound annotations]',
                    help='path to genome bound file to limit the chromosome for analysis')
    args = parser.parse_args()
    # input argument variables
    GENCODE = args.gencode_annotations
    out_gene_bed = args.out_gene_bed
    out_exon_bed = args.out_exon_bed
    out_gene_tss = args.out_gene_tss
    out_gene_tes = args.out_gene_tes
    genome_bound = args.genome_bound
    # process files. 
    GENCODE_df = pr.read_gtf(GENCODE)
    GENCODE_df=GENCODE_df.as_df()
    genome_bound_df = pd.read_csv(genome_bound,sep='\t',header=None,names=['chrom','length'])
    chromosomes = genome_bound_df['chrom']
    # filter to only chromosomes in bound file
    GENCODE_df=GENCODE_df[GENCODE_df['Chromosome'].isin(chromosomes)]
    # genes
    condition_v7_gene_range = (GENCODE_df['Feature']=='gene')
    v7_gene_range_df = GENCODE_df[condition_v7_gene_range].copy()
    v7_gene_range_df['Bed_start'] = v7_gene_range_df['Start']-1
    genes=v7_gene_range_df[(v7_gene_range_df['gene_type']=='protein_coding')|(v7_gene_range_df['gene_type']=='lincRNA')][['Chromosome','Bed_start','End','gene_id','gene_type','Strand']]

    genes[~genes.gene_id.str.contains('PAR')].to_csv(out_gene_bed,header=False,index=False,sep='\t')
    # exons
    condition_v7_exon_range = (GENCODE_df['Feature']=='exon')
    v7_exon_range_df = GENCODE_df[condition_v7_exon_range].copy()
    v7_exon_range_df['Bed_start'] = v7_exon_range_df['Start']-1
    exons=v7_exon_range_df[(v7_exon_range_df['gene_type']=='protein_coding')|(v7_exon_range_df['gene_type']=='lincRNA')][['Chromosome','Bed_start','End','gene_id','exon_number']]
    exon_pl=pl.from_pandas(exons)
    exon_pl=exon_pl.with_columns([(pl.col('End')-pl.col('Bed_start')).alias('exon_length'),pl.col('exon_number').count().over('gene_id').alias('num_exons')])
    exon_pl=exon_pl.with_columns(pl.col('exon_length').sum().over('gene_id').alias('total_coding_length'))
    exon_pl.write_csv(out_exon_bed,has_header=False,separator='\t')
    # tss tes
    genes['tss_start'] = genes.apply(lambda x: x.Bed_start if (x.Strand=='+') else x.End-1,axis=1)
    genes['tss_end'] = genes.apply(lambda x: x.Bed_start+1 if (x.Strand=='+') else x.End,axis=1)
    genes['tes_start'] = genes.apply(lambda x: x.End-1 if (x.Strand=='+') else x.Bed_start,axis=1)
    genes['tes_end'] = genes.apply(lambda x: x.End if (x.Strand=='+') else x.Bed_start+1,axis=1)

    genes[['Chromosome','tss_start','tss_end','gene_id']].to_csv(out_gene_tss,header=False,index=False,sep='\t')
    genes[['Chromosome','tes_start','tes_end','gene_id']].to_csv(out_gene_tes,header=False,index=False,sep='\t')

