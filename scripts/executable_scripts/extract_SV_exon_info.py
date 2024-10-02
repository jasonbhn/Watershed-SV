#!/usr/bin/env python3

import pandas as pd
import polars as pl
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gene to SV distance')
    parser.add_argument('--input',type=str,required=True,metavar='[input file path]',
                    help='input file path')
    parser.add_argument('--output',type=str,required=True,metavar='[output file path]',
                    help='output file path')
    args = parser.parse_args()
    # input argument variables
    inputfile = args.input
    outputfile = args.output
    exon_info_unprocessed = pl.scan_csv(inputfile,separator='\t',
                has_header=False,
                new_columns=['Gene','SV','exon','exon_length','total_exon_num','total_exon_length','exon_overlap'],
                dtypes={'SV':str})
    exon_info_unprocessed = exon_info_unprocessed.with_columns(
        (pl.col('exon_overlap')/pl.col('exon_length')).alias('exon_fraction_affected'))
    exon_info_processed = exon_info_unprocessed.groupby(['Gene','SV']).agg([
        (pl.col('exon_overlap').sum()/pl.col('total_exon_length')).first().alias('coding_fraction_affected'),
        pl.when(pl.col('exon_fraction_affected').first()<1).then(1).otherwise(0).alias('SV_5prime_exon_truncation'),
        pl.when(pl.col('exon_fraction_affected').last()<1).then(1).otherwise(0).alias('SV_3prime_exon_truncation')
    ]).collect()
    exon_info_processed = exon_info_processed.with_columns(
        pl.when(pl.col('coding_fraction_affected')!=0).then(1).otherwise(0).alias('exon_variant')
    )
    exon_info_processed.write_csv(outputfile,separator='\t',has_header=True,null_value='')
