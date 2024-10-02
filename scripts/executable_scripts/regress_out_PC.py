#!/usr/bin/env python3

import os
import glob
import mygene
import numpy as np
import pandas as pd
import seaborn as sns
import pysam
import gseapy
import argparse
import lxml
from pca import pca
from scipy import stats
from sklearn.datasets import load_iris
from matplotlib import pyplot as plt
from statsmodels.api import OLS
from dask import delayed
from dask import dataframe
from dask.distributed import Client,LocalCluster
def unit_regressor_disease_specific(gene,y_data,expression_PC):
    y = y_data
    X = expression_PC
    df=pd.concat([y,X],axis=1)
    expression_pcs = [i for i in expression_PC.columns if 'PC' in i]
    covars = [i for i in expression_PC.columns if 'PC' not in i]
    model = OLS.from_formula(formula=f'{gene} ~ {" + ".join(expression_pcs)}',data=df)
    result = model.fit()
    category_to_param_key = {i:f'C(Disease_class)[T.{i}]' for i in X['Disease_class'].unique().tolist()}
    unnorm_corrected = list(result.resid+[result.params[category_to_param_key[i]] if i in result.params else 0 for i in X['Disease_class']])
    corrected = (unnorm_corrected-np.mean(unnorm_corrected))/np.std(unnorm_corrected)
    return pd.DataFrame({'gene':gene,
                         'Ind':y.index,
                         'normalized_resid':corrected}).astype({'gene':'string','Ind':'string','normalized_resid':'float64'})

def unit_regressor(gene,y_data,expression_PC):    
    y = y_data
    X = expression_PC
    df=pd.concat([y,X],axis=1)
    expression_pcs = [i for i in expression_PC.columns if 'PC' in i]
    covars = [i for i in expression_PC.columns if 'PC' not in i]
    model = OLS.from_formula(formula=f'{gene} ~ {" + ".join(expression_pcs)} + C(Sex)',data=df)
    result = model.fit()
    unnorm_corrected = list(result.resid)
    corrected = (unnorm_corrected-np.mean(unnorm_corrected))/np.std(unnorm_corrected)
    return pd.DataFrame({'gene':gene,
                         'Ind':y.index,
                         'normalized_resid':corrected}).astype({'gene':'string','Ind':'string','normalized_resid':'float64'})

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='VCF filtering based on allele frequency field')
    parser.add_argument('--regression-data',type=str,required=True,metavar='[input regressiond data file]',
                    help='select data.txt to parse')
    parser.add_argument('--num-PC-remove',type=int,metavar='[remove n number of top PCs]',
                    default=1,
                    help='how many top PCs to remove from the data')
    parser.add_argument('--out-parquet-dir',type=str,required=True,metavar='[input for downstream annotSV]',
                    help='path to input file to be created for annotSV')
    args = parser.parse_args()
    regression_data_w_sex=pd.read_csv(args.regression_data,sep='\t',index_col=0)

    template_meta = pd.DataFrame({'gene': pd.Series(dtype='string'),
                              'Ind': pd.Series(dtype='string'),
                              'normalized_resid':pd.Series(dtype='float64')})
    residuals=[]
    gene_cols = regression_data_w_sex.columns[regression_data_w_sex.columns.str.contains('ENSG')]
    hidden_confound_cols = regression_data_w_sex.columns[regression_data_w_sex.columns.str.contains('PC')].tolist()
    PCs_to_regress_cols = hidden_confound_cols[:args.num_PC_remove]
    known_covars_cols = regression_data_w_sex.columns[(~regression_data_w_sex.columns.str.contains('PC'))&(~regression_data_w_sex.columns.str.contains('ENSG'))].tolist()
    for gene in gene_cols:
        delaycompute=delayed(unit_regressor_disease_specific)(gene,
                                            regression_data_w_sex[gene],
                                            regression_data_w_sex[PCs_to_regress_cols+known_covars_cols])
        residuals.append(delaycompute)
    collector=dataframe.from_delayed(residuals,meta=template_meta)
    repartitioned_collector=collector.repartition(npartitions=6)
    repartitioned_collector.to_parquet(args.out_parquet_dir,compression='snappy')
