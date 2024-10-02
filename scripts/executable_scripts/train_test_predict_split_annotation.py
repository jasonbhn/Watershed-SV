#!/usr/bin/env python3

from pysam import VariantFile
from sv_utils import *
import pandas as pd
import functools
import os
import numpy as np
import argparse

def drop_uninformative_columns(dataframe,mode='evaluate'):
    to_drop=[]
    if mode=='evaluate':
        training=dataframe[dataframe.N2pair=="NA"]
    elif mode=='predict':
        training=dataframe
    for i in training.columns:
        if i!="N2pair" and i!="partition" and training[i].value_counts().max()==training.shape[0]:
            to_drop.append(i)
            print(i)
    return dataframe.drop(to_drop,axis=1),to_drop

# again this method has many modes:
# mode 1: split a data into n2-unique set, which is the input
# for evaluate_watershed.R
# mode 2: given two or more datasets, where one is intended as the 
# training data, the others as the test data, jointly QC and clean

if __name__ == '__main__':
    # parse the arguments. 
    parser = argparse.ArgumentParser(description='dataset utility')
    parser.add_argument('--training',type=str,required=True,metavar='[input training data]',
                    help='which training dataset')
    parser.add_argument('--testings',nargs='+',default=[],metavar='[input one or more test data]',
                    help='which are the testing data')
    parser.add_argument('--testing-list',type=str,default='None',metavar='[txt contain files to merge for test]',
                    help='a txt file each line is path to a test file. ')
    parser.add_argument('--gene-ind-sv',type=str,default='None',metavar='[SV list per gene-ind pair]',
                    help='file to help select N2 pairs')
    parser.add_argument('--mode',required=True,choices=['evaluate', 'predict'],metavar='[choose how data should be generated]',
                    help='options corresponds to either generating output for predict or evaluate Watershed')
    parser.add_argument('--min-af-impute-mode',required=True,choices=['infer', 'upload'],metavar='[choose how to impute af=0]',
                    help='options corresponds to either picking minimum af in data or manually enter impute value')
    parser.add_argument('--min-af-value',default=1/831,type=float,metavar='[value to impute for af=0]',
                    help='if option is upload then use this value.')
    parser.add_argument('--out-prefix',type=str,required=True,metavar='[output cleaned data]',
                    help='prefix, including file path to prefix')
    args = parser.parse_args()

    if args.mode=='predict':
        if args.testing_list=='None' and len(args.testings)==0:
            print('mode==predict, please provide one of --testings or --testing-list')
            exit()
        training=pd.read_csv(args.training,sep='\t')
        test_dfs = []
        if args.testing_list != 'None':
            with open(args.testing_list,'r') as handle:
                for i in handle.readlines():
                    test_dfs.append(pd.read_csv(i.rstrip(),sep='\t'))
        else:
            for i in args.testings:
                test_dfs.append(pd.read_csv(i,sep='\t'))
        testing=pd.concat(test_dfs)
        if 'SVTYPE_INV' in testing:
            testing.loc[testing['SVTYPE_INV']==1,'SVTYPE_DEL'] = 1
            testing.loc[testing['SVTYPE_INV']==1,'SVTYPE_INS'] = 1
        if 'SVTYPE_MEI' in testing:
            testing.loc[testing['SVTYPE_MEI']==1,'SVTYPE_DEL'] = 1
        # add partition marker
        training['partition'] = 'train'
        testing['partition'] = 'test'
        to_drop = list(set(testing.columns).difference(set(training.columns)))
        total_data = pd.concat([training,testing]).drop(to_drop,axis=1)
        outlier_cols = total_data.columns[total_data.columns.str.contains('TE_pvalues')].to_list()
        total_data['N2pair'] = 'NA'
        total_data[outlier_cols] = total_data[outlier_cols].fillna('NaN')
        total_data.af.fillna(0,inplace=True)
        # impute af=0
        if args.min_af_impute_mode == 'upload':
            total_data['af'] = np.where(total_data['af']==0,args.min_af_value,total_data['af'])
        else:
            inferred_min = total_data[total_data['af']!=0].af.min()
            total_data['af'] = np.where(total_data['af']==0,inferred_min,total_data['af'])
        # impute all the rest
        methods=pd.read_csv('/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/collapse_annotation_instructions.tsv',sep='\t')
        impute_records=methods[~methods.Impute_method.isna()].to_records()
        impute_dict={record[1]:record[2] for record in impute_records}
        for column in total_data.columns:
            if column in impute_dict:
                if impute_dict[column]=='0':
                    total_data[column]=total_data[column].fillna(0)
                elif impute_dict[column]=='mean':
                    impute_val = total_data.loc[total_data.N2pair=='NA',column].mean()
                    total_data[column]=total_data[column].fillna(impute_val)
        total_data.fillna(0,inplace=True)
        # then drop the uninformative columns from the data. 
        _,drop_cols = drop_uninformative_columns(total_data[total_data.partition=='train'],mode='predict')
        total_data[total_data.partition=="train"][[c for c in total_data if c not in drop_cols +['partition', 'N2pair']+outlier_cols] + outlier_cols + ['N2pair']].\
        to_csv(f'{args.out_prefix}_training_data.tsv',sep='\t',header=True,index=False)
        total_data[total_data.partition=="test"][[c for c in total_data if c not in drop_cols+['partition','N2pair']+outlier_cols] + outlier_cols + ['N2pair']].\
        to_csv(f'{args.out_prefix}_testing_data.tsv',sep='\t',header=True,index=False)


                
