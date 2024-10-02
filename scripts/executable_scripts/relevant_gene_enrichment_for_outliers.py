#!/usr/bin/env python3

import argparse
import os
import glob
import numpy as np
import pandas as pd
import polars as pl
import pysam
from scipy import stats

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='script to generate rare SV enrichment')
    parser.add_argument('--input-file',type=str,required=True,metavar='[input residuals file paths]',
                    help='residual file to process')
    parser.add_argument('--relevant-gene-set-file',type=str,required=True,metavar='[input relevant gene set file paths]',
                    help='relevant gene set to process. ')
    parser.add_argument('--num-PCs',type=int,required=True,metavar='[number of PCs removed]',
                    help='number of PCs removed')
    parser.add_argument('--outfile',type=str,required=True,metavar='[output the enrichment statistics]',
                    help='output enrichment stat. ')
    args = parser.parse_args()
    filename = os.path.join(args.input_file,'*.parquet')
    out = args.outfile
    num_PCs = args.num_PCs
    relevant_gene_set_file = args.relevant_gene_set_file
    # read in list of ensembl ids for relevant gene set
    opentarget_ensembl_list = pd.read_csv(relevant_gene_set_file,sep='\t',header=None,names=['gene'])['gene'].tolist()
    # iterate through all PCs to calculate rare SV enrichment
    pre_outliers = pl.scan_parquet(filename)
    pre_outliers = pre_outliers.with_columns([pl.when(pl.col('normalized_resid').abs() > 3).then(1).otherwise(0).alias('Y'),
                                              pl.when(pl.col('Ind').str.contains('GTEX')).then('GTEx').otherwise('CMG').alias('cohort'),
                                              pl.col('gene').str.split('.').arr[0]
                                             ])

    relevant_outliers = pre_outliers.with_columns(pl.when(pl.col('gene').is_in(opentarget_ensembl_list)).then(True).otherwise(False).alias('MuscleDisease_relevant'))
    relevant_outliers = relevant_outliers.with_columns([pl.when(pl.col('normalized_resid').abs() > i).then(1).otherwise(0).alias(f"Y{i}") for i in range(1,6)]).collect()


    cmg_enriched_relevant_genes = relevant_outliers.filter(pl.col('cohort')=='CMG').\
    select([pl.concat_list(pl.col([f'Y{i}','MuscleDisease_relevant'])).alias('Y_relevant_category').
            value_counts(sort=True).struct.rename_fields(['Y_relevant_category','count']).
            alias(f'Gene_cat_count_{i}') for i in range(1,6)])
    cmg_contingency_tables=[cmg_enriched_relevant_genes.select(pl.col(f'Gene_cat_count_{i}').struct.field('Y_relevant_category').cast(pl.List(str)).arr.join(","),
                                pl.col(f'Gene_cat_count_{i}').struct.field('count').alias(f'count_Y{i}'),
                                )for i in range(1,6)]
    cmg_contingency_table = cmg_contingency_tables[0]
    for i in range(1,5):
        cmg_contingency_table = cmg_contingency_table.join(cmg_contingency_tables[i],on='Y_relevant_category')
    cmg_contingency_table_transpose=cmg_contingency_table.sort(by="Y_relevant_category").select(pl.all().exclude('Y_relevant_category')).\
    transpose(column_names=['0,0','0,1','1,0','1,1'])
    cmg_relevant_genes_enrichment = cmg_contingency_table_transpose.with_columns([pl.Series(values=[1,2,3,4,5]).alias('Z-threshold'),
                                            (pl.col('1,1')/pl.col('0,1')/((pl.col('1,0')/pl.col('0,0')))).alias('enrichment'),
                                            pl.struct(['1,1','0,1','1,0','0,0']).
                                            apply(lambda x: stats.fisher_exact(np.array([[x['1,1'],x['0,1']],
                                                                                        [x['1,0'],x['0,0']]]),
                                                                                alternative='two-sided')[1]).alias('p-value'),
                                                                                    pl.lit('CMG').alias("cohort"),
                                                                                pl.lit(num_PCs).alias("number_of_PCs")
                                            ])
    gtex_enriched_relevant_genes = relevant_outliers.filter(pl.col('cohort')=='GTEx').\
    select([pl.concat_list(pl.col([f'Y{i}','MuscleDisease_relevant'])).alias('Y_relevant_category').
            value_counts(sort=True).struct.rename_fields(['Y_relevant_category','count']).
            alias(f'Gene_cat_count_{i}') for i in range(1,6)])
    gtex_contingency_tables=[gtex_enriched_relevant_genes.select(pl.col(f'Gene_cat_count_{i}').struct.field('Y_relevant_category').cast(pl.List(str)).arr.join(","),
                                pl.col(f'Gene_cat_count_{i}').struct.field('count').alias(f'count_Y{i}'),
                                )for i in range(1,6)]
    gtex_contingency_table = gtex_contingency_tables[0]
    for i in range(1,5):
        gtex_contingency_table = gtex_contingency_table.join(gtex_contingency_tables[i],on='Y_relevant_category')
    gtex_contingency_table_transpose=gtex_contingency_table.sort(by="Y_relevant_category").select(pl.all().exclude('Y_relevant_category')).\
    transpose(column_names=['0,0','0,1','1,0','1,1'])
    gtex_relevant_genes_enrichment = gtex_contingency_table_transpose.with_columns([pl.Series(values=[1,2,3,4,5]).alias('Z-threshold'),
                                            (pl.col('1,1')/pl.col('0,1')/((pl.col('1,0')/pl.col('0,0')))).alias('enrichment'),
                                            pl.struct(['1,1','0,1','1,0','0,0']).
                                            apply(lambda x: stats.fisher_exact(np.array([[x['1,1'],x['0,1']],
                                                                                        [x['1,0'],x['0,0']]]),
                                                                                alternative='two-sided')[1]).alias('p-value'),
                                                                                    pl.lit('GTEx').alias("cohort"),
                                                                                pl.lit(num_PCs+1).alias("number_of_PCs")
                                            ])
    pl.concat([cmg_relevant_genes_enrichment,gtex_relevant_genes_enrichment]).write_csv(out,separator='\t',has_header=True)
    print(f'calculated_enrichment for removing {num_PCs} PCs')
