import sys
new_path = '/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/scripts/executable_scripts/'

if new_path not in sys.path:
    sys.path.append(new_path)
from pysam import VariantFile
from sv_utils import *
import pandas as pd
import pyranges as pr
import polars as pl
import functools
import shutil
import numpy as np
import os
import argparse
import polars.selectors as cs 

def drop_uninformative_columns(dataframe,mode='evaluate'):
    to_drop=[]
    if mode=='evaluate':
        training=dataframe[dataframe.N2pair=="NA"]
    elif mode=='predict':
        training=dataframe
    for i in training.columns:
        if i!="N2pair"and training[i].value_counts().max()==training.shape[0]:
            to_drop.append(i)
            print(i)
    return dataframe.drop(to_drop,axis=1),to_drop

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='tools to clean up, impute watershed input')
    parser.add_argument('--gene-sv-annotation',type=str,required=True,metavar='[gene-sv annot for n2]',
                    help='gene sv annotations')
    parser.add_argument('--gene-annotation',type=str,required=True,metavar='[gene annot for n2]',
                    help='gene annotations')
    parser.add_argument('--collapse-instructions',type=str,required=True,metavar='[how to collapse and impute]',
                    help='gene annotations imputation instructions')
    parser.add_argument('--pval-threshold',type=float,required=False,default=0.0027, metavar='[p-value threshold]',
                    help='select p-value threshold')
    parser.add_argument('--output',type=str,required=True,metavar='[output eval frame]',
                    help='output eval frame')
    args = parser.parse_args()

    # io some file. 
    gene_sv_annotation = pd.read_csv(args.gene_sv_annotation,sep='\t')
    gene_annotation = pd.read_csv(args.gene_annotation,sep='\t')
    methods=pd.read_csv(args.collapse_instructions,sep='\t')
    pthreshold = args.pval_threshold
    gene_sv_annotation[['GeneName','SV']] = gene_sv_annotation['GeneName'].str.split(':',expand=True)
    # get n2 pair samples. 
    for_n2=gene_sv_annotation[['SubjectID','GeneName','SV','dCN','Allele']]
    for_n2['SV_dCN_Alleles'] = for_n2.SV + ':' + for_n2.dCN.astype(str) + ':' + for_n2.Allele.astype(str)
    grouping=for_n2.groupby(['SubjectID','GeneName']).agg(SVs=('SV_dCN_Alleles',lambda x: ','.join(np.sort(x.unique()))),
                                                        dCN=('dCN',lambda x: x.loc[x.abs().idxmax()]),
                                                        Allele=('Allele','max')).reset_index()
    out=[]
    n2_id=0
    predict_out=[]
    training_out=[]

    for index,df in grouping.groupby(['SVs','GeneName']):
        if df.shape[0]>1:
            n2pair=df.sample(n=2)
            n2pair['N2pair']=n2_id
            df['N2pair']=n2_id
            n2_id+=1
            out.append(n2pair)
            predict_out.append(df)
        else:
            df['N2pair']='NA'
            out.append(df)
            training_out.append(df)

    training_testing=pd.concat(out)
    predict_df=pd.concat(predict_out).drop('N2pair',axis=1)
    training_df=pd.concat(training_out).drop('N2pair',axis=1)
    train_test_split=training_testing.drop(['SVs'],axis=1)
    predict_df.drop(['SVs'],axis=1,inplace=True)
    training_df.drop(['SVs'],axis=1,inplace=True)
    # fetch list of outlier columns. 
    pval_list = list(gene_annotation.columns[gene_annotation.columns.str.contains('TE_pvalues')])
    # split to train test then recombine.
    output_benchmark_train = gene_annotation.merge(train_test_split.loc[train_test_split['N2pair']=='NA',['SubjectID','GeneName','N2pair']],
                        on=['SubjectID','GeneName'])[[i for i in gene_annotation.columns if 'TE_pvalues' not in i ]+pval_list+['N2pair']]

    output_benchmark_test = gene_annotation.merge(train_test_split.loc[train_test_split['N2pair']!='NA',['SubjectID','GeneName','N2pair']],
                        on=['SubjectID','GeneName'])[[i for i in gene_annotation.columns if 'TE_pvalues' not in i ]+pval_list+['N2pair']]
    output_benchmark_total = pd.concat([output_benchmark_train,output_benchmark_test])

    # print the proportion outlier value
    outlier_fracs = []
    for p in pval_list:
        frac = (output_benchmark_total[p].abs()<pthreshold).value_counts(normalize=True)
        outlier_fracs.append(frac)
        print(f"outlier fraction, {p}: {frac}")

    # imputations
    impute_records=methods[~methods.Impute_method.isna()].to_records()
    impute_dict={record[1]:record[2] for record in impute_records}
    collapse_records=methods.to_records()
    collapse_dict={record[1]:record[3] for record in collapse_records}
    del collapse_dict['GeneName']
    del collapse_dict['SubjectID']

    # get pvalues:
    collapse_modified = {}
    for  key,val in collapse_dict.items():
        if key in output_benchmark_total.columns and key !='dCN':
            collapse_modified[key]=val
        elif key == 'dCN':
            collapse_modified[key]=(lambda x:x.loc[x.abs().idxmax()])

    
    output_benchmark_total[pval_list] = output_benchmark_total[pval_list].fillna('NaN')

    output_benchmark_total[output_benchmark_total.columns[~output_benchmark_total.columns.isin(pval_list)]] = output_benchmark_total[output_benchmark_total.columns[~output_benchmark_total.columns.isin(pval_list)]].fillna(0)

    cleaned_benchmark, columns_dropped = drop_uninformative_columns(output_benchmark_total)

    pd.concat([cleaned_benchmark[cleaned_benchmark.N2pair=='NA'],
           cleaned_benchmark[cleaned_benchmark.N2pair!='NA'].sort_values(by='N2pair')]).\
    to_csv(args.output,
            sep='\t',
            header=True,
            index=False)
