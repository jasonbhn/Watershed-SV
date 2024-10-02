#!/usr/bin/env python3

import polars as pl
import polars.selectors as cs
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='curate SV specific annotations')
    parser.add_argument('--vcf',type=str,required=True,metavar='[input rare SV VCF file]',
                    help='select [rare_SV].vcf to parse')
    parser.add_argument('--genotypes',type=str,required=True,metavar='[input genotype VCF file]',
                    help='select [rare_SV].genotype to parse')
    parser.add_argument('--genes',type=str,required=True,metavar='[input gene list file]',
                    help='select genes to parse')

    gene_sv_gb_df = pl.scan_csv(gene_sv_gb,
                            separator='\t',
                            has_header=False,
                            new_columns=['chrom','start','end','SV','SVTYPE','Gene'])
    gene_sv_tss_df = pl.scan_csv(gene_sv_tss,
                            separator='\t',
                            has_header=False,
                            new_columns=['chrom','start','end','SV','SVTYPE','Gene'])
    gene_sv_tes_df = pl.scan_csv(gene_sv_tes,
                            separator='\t',
                            has_header=False,
                            new_columns=['chrom','start','end','SV','SVTYPE','Gene'])
