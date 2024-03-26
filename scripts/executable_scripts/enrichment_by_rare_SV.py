import argparse
import os
import glob
import mygene
import numpy as np
import pandas as pd
import polars as pl
import seaborn as sns
import pysam
from scipy import stats
from matplotlib import pyplot as plt
from statsmodels.api import OLS

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='script to generate rare SV enrichment')
    parser.add_argument('--input-file-list',type=str,required=True,metavar='[input residuals file paths]',
                    help='list of residuals to process in a txt, in order of number of PCs')
    parser.add_argument('--outfile',type=str,required=True,metavar='[output the enrichment statistics]',
                    help='output of the enrichment analysis')
    parser.add_argument('--sv-annot-dir',type=str,required=True,metavar='[input genotype file paths]',
                    help='directory containing all watershed-sv annotations')
    parser.add_argument('--maf-file',type=str,required=True,metavar='[input maf file paths]',
                    help='maf file path')
    parser.add_argument('--flank',type=int,required=True,metavar='[flanking distance]',
                    help='flanking distance')
    parser.add_argument('--af-colname',type=str,required=True,metavar='[colname for af]',
                    help='colname for af')
    parser.add_argument('--sample-colname',type=str,required=True,metavar='[colname for sample]',
                    help='colname for sample')
    parser.add_argument('--expression-colname',type=str,required=True,metavar='[colname for expression]',
                    help='colname for expression')
    parser.add_argument('--use-gtex',action=argparse.BooleanOptionalAction,metavar='[does expression data combine with gtex?]',
                    help='use gtex or not?')
    args = parser.parse_args()
    maf_file = args.maf_file
    af_colname = args.af_colname
    sample_colname = args.sample_colname
    expression_colname = args.expression_colname
    use_gtex = args.use_gtex
    sv_annot_dir = args.sv_annot_dir
    flank = args.flank
    with open(args.input_file_list,'r') as handle:
        files_newline = handle.readlines()
    if 'pq'in files_newline[0] or 'parquet' in files_newline[0]:
        file_type='parquet'
        files = [i.rstrip() for i in files_newline]
    elif 'txt' in files_newline[0]:
        file_type='txt'
        files = [i.rstrip() for i in files_newline]
    else:
        file_type='parquet'
        files = [os.path.join(i.rstrip(),'*.parquet') for i in files_newline]
    # iterate through all PCs to calculate rare SV enrichment
    total_enrichment_by_number_of_PCs = []
    count = 0
    for expression_file in files:
        if file_type =='parquet':
            pre_outliers = pl.scan_parquet(expression_file)
        else:
            pre_outliers = pl.scan_csv(expression_file,separator='\t')
        pre_outliers = pre_outliers.with_columns([pl.when(pl.col(expression_colname).abs() > 3).then(1).otherwise(0).alias('Y'),
                                                pl.when(pl.col(sample_colname).str.contains('GTEX')).then('GTEx').otherwise('CMG').alias('cohort'),
                                                pl.col('gene').str.split('.').arr[0]
                                                ])
        SV_genotypes_cmg=pl.scan_csv(os.path.join(sv_annot_dir,'intermediates/pipeline_input_genotypes.tsv'),separator='\t',dtypes={'SV':str})
        gene_SV_overlap_cmg=pl.scan_csv(os.path.join(sv_annot_dir,f'intermediates/gene_sv.{flank}.bed'),separator='\t',dtypes={'SV':str},
                                new_columns=['chrom','start','end','SV','SVTYPE','gene'])
        gene_SV_overlap_cmg=gene_SV_overlap_cmg.with_columns(pl.col('gene').str.split('.').arr.first())
        if use_gtex:
            SV_genotypes_gtex=pl.scan_csv('new-protein-lincRNA-10k-8.0/intermediates/pipeline_input_genotypes.tsv',separator='\t',dtypes={'SV':str})
            gene_SV_overlap_gtex=pl.scan_csv('new-protein-lincRNA-10k-8.0/intermediates/gene_sv.10000.bed',separator='\t',dtypes={'SV':str},
                                new_columns=['chrom','start','end','SV','SVTYPE','gene'])
            gene_SV_overlap_gtex=gene_SV_overlap_gtex.with_columns(pl.col('gene').str.split('.').arr.first())
            GTEX_MAF = pl.scan_csv('new-protein-lincRNA-10k-8.0/intermediates/custom_gtex_maf.tsv',separator='\t',)
                    
            rare_svs_gtex=GTEX_MAF.filter(pl.col('af')<0.01).select('SV').collect().to_series()
            rare_genotypes_gtex = SV_genotypes_gtex.filter((pl.col('Allele')>0)&
                                (pl.col('SV').is_in(rare_svs_gtex))).collect().rename({'SUBJID':'Ind'}).drop('SVid').lazy()
        CMG_MAF = pl.scan_csv(maf_file,separator='\t',dtypes={'SV':str})
        rare_svs_cmg=CMG_MAF.filter(pl.col(af_colname)<0.01).select('SV').collect().to_series()
        rare_genotypes_cmg = SV_genotypes_cmg.filter((pl.col('Allele')>0)&
                            (pl.col('SV').is_in(rare_svs_cmg))).collect().rename({'SUBJID':sample_colname}).lazy()

        relevant_outliers=pre_outliers.with_columns([pl.when(pl.col(expression_colname).abs() > i).then(1).otherwise(0).alias(f"Y{i}") for i in range(1,6)])
        rare_sv_enrichment_cmg = relevant_outliers.join(gene_SV_overlap_cmg.select(['gene','SV']),on='gene',how='left').\
        join(rare_genotypes_cmg,on=[sample_colname,'SV'],how='left').filter(pl.col('cohort')=='CMG').\
        with_columns(pl.when(pl.col('Allele').is_null()).then(0).otherwise(1).alias('has_rare_SV')).collect()
        if use_gtex:
            # both cmg and gtex together. 
            gene_SV_overlap_both = pl.concat([gene_SV_overlap_cmg,gene_SV_overlap_gtex])
            rare_genotypes_both = pl.concat([rare_genotypes_cmg,rare_genotypes_gtex])
            rare_sv_enrichment_both = relevant_outliers.join(gene_SV_overlap_both.select(['gene','SV']),on='gene',how='left').\
            join(rare_genotypes_both,on=['Ind','SV'],how='left').\
            with_columns(pl.when(pl.col('Allele').is_null()).then(0).otherwise(1).alias('has_rare_SV')).collect()

            both_enriched_sv_genes = rare_sv_enrichment_both.\
            select([pl.concat_list(pl.col([f'Y{i}','has_rare_SV'])).alias('Y_rare_SV').
                    value_counts(sort=False).struct.rename_fields(['Y_rare_SV','count']).
                    alias(f'Gene_cat_count_{i}') for i in range(1,6)])
        
            both_contingency_tables=[both_enriched_sv_genes.select(pl.col(f'Gene_cat_count_{i}').struct.field('Y_rare_SV').cast(pl.List(str)).arr.join(","),
                                    pl.col(f'Gene_cat_count_{i}').struct.field('count').alias(f'count_Y{i}'),
                                    )for i in range(1,6)]

            both_contingency_table = both_contingency_tables[0]
            for i in range(1,5):
                both_contingency_table = both_contingency_table.join(both_contingency_tables[i],on='Y_rare_SV')

            both_contingency_table_transpose=both_contingency_table.sort(by="Y_rare_SV").select(pl.all().exclude('Y_rare_SV')).\
            transpose(column_names=['0,0','0,1','1,0','1,1'])

            both_rare_sv_genes_enrichment = both_contingency_table_transpose.with_columns([pl.Series(values=[1,2,3,4,5]).alias('Z-threshold'),
                                                (pl.col('1,1')/pl.col('0,1')/((pl.col('1,0')/pl.col('0,0')))).alias('enrichment'),
                                                pl.struct(['1,1','0,1','1,0','0,0']).
                                                apply(lambda x: stats.fisher_exact(np.array([[x['1,1'],x['0,1']],
                                                                                            [x['1,0'],x['0,0']]]),
                                                                                    alternative='greater')[1]).alias('p-value'),
                                                                                        pl.lit('GTEx+CMG').alias("cohort"),
                                                                                    pl.lit(count+1).alias("number_of_PCs")
                                                ])
        cmg_enriched_sv_genes = rare_sv_enrichment_cmg.filter(pl.col('cohort')=='CMG').\
        select([pl.concat_list(pl.col([f'Y{i}','has_rare_SV'])).alias('Y_rare_SV').
                value_counts(sort=False).struct.rename_fields(['Y_rare_SV','count']).
                alias(f'Gene_cat_count_{i}') for i in range(1,6)])

        cmg_contingency_tables=[cmg_enriched_sv_genes.select(pl.col(f'Gene_cat_count_{i}').struct.field('Y_rare_SV').cast(pl.List(str)).arr.join(","),
                                    pl.col(f'Gene_cat_count_{i}').struct.field('count').alias(f'count_Y{i}'),
                                    )for i in range(1,6)]

        cmg_contingency_table = cmg_contingency_tables[0]
        for i in range(1,5):
            cmg_contingency_table = cmg_contingency_table.join(cmg_contingency_tables[i],on='Y_rare_SV')
        cmg_contingency_table_transpose=cmg_contingency_table.sort(by="Y_rare_SV").select(pl.all().exclude('Y_rare_SV')).\
        transpose(column_names=['0,0','0,1','1,0','1,1'])

        cmg_rare_sv_genes_enrichment = cmg_contingency_table_transpose.with_columns([pl.Series(values=[1,2,3,4,5]).alias('Z-threshold'),
                                                (pl.col('1,1')/pl.col('0,1')/((pl.col('1,0')/pl.col('0,0')))).alias('enrichment'),
                                                pl.struct(['1,1','0,1','1,0','0,0']).
                                                apply(lambda x: stats.fisher_exact(np.array([[x['1,1'],x['0,1']],
                                                                                            [x['1,0'],x['0,0']]]),
                                                                                    alternative='greater')[1]).alias('p-value'),
                                                                                        pl.lit('CMG').alias("cohort"),
                                                                                    pl.lit(count+1).alias("number_of_PCs")
                                                ])
        if use_gtex:
            total_enrichment_by_number_of_PCs.append(pl.concat([cmg_rare_sv_genes_enrichment,both_rare_sv_genes_enrichment]))
        else:
            total_enrichment_by_number_of_PCs.append(cmg_rare_sv_genes_enrichment)
        print(f'calculated_enrichment for removing {count+1} PCs')
        count+=1
    pl.concat(total_enrichment_by_number_of_PCs).\
    write_parquet(args.outfile)
