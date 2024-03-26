import numpy as np
import pandas as pd
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
import scipy.stats as stats
import statsmodels as sm
from statsmodels.stats.multitest import fdrcorrection
from sklearn.impute import KNNImputer
import seaborn as sns
import pickle
from io import StringIO
import glob
import functools
from matplotlib import pyplot as plt

def cn_change(svtype):
    if svtype=='DEL_CNV':
        return '-'
    elif svtype=='DUP_CNV':
        return '+'
    elif svtype=='DUP':
        return '+'
    elif svtype=='DEL':
        return '-'
    elif svtype=='BND':
        return '='
    elif svtype=='INV':
        return '='
    elif svtype=='ALU':
        return '+'
    elif svtype=='LINE1':
        return '+'
    elif svtype=='MEI':
        return '-'
    elif svtype=='SVA':
        return '+'

def serialize_indiv_sv(ac_matrix):
    cat = []
    for i in ac_matrix.columns[:-1]:
        SVs = ac_matrix[i][ac_matrix[i]!=0].index.tolist()
        cat.extend(list(tuple(zip([i]*len(SVs),SVs))))
    return cat

def Zscores_to_Pvalues(Zscores):
    return stats.norm.sf(abs(Zscores))*2*np.sign(Zscores)

gene_indiv_set=pd.read_csv('starting_annotations/medZ_gene_indiv_set.csv',sep='\t')
annotations=pd.read_csv('tissue_specific_raw_annotations.10000.bed',sep='\t')
remap_addition=pd.read_csv('starting_annotations/overlaps_remap_crm.10000.tsv',sep='\t')
var_matrix=pd.read_csv('variants/gtex.lumpy.gs.melt.high_conf.hg38.acMatrix.tsv',sep='\t')
sample_phenotypes=pd.read_csv('/work-zfs/abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt',sep='\t')
var_matrix=var_matrix.set_index('Unnamed: 0')

annotations.drop(['chrom','start','end'],axis=1,inplace=True)
annotations[[i for i in annotations.columns if 'Z:' not in i]]=annotations[[i for i in annotations.columns if 'Z:' not in i]].fillna(0)
avail_tissues = annotations[[i for i in annotations.columns if 'Z:' in i]].count(axis=1)
annot_5_tissues = annotations[avail_tissues >=5 ] # july-changes to better match multitissue outlier signals. 
annot_5_tissues['Y']=annot_5_tissues[[i for i in annot_5_tissues.columns if 'Z:' in i]].apply(lambda x: 'outlier' if (abs(np.nanmedian(x.to_numpy())) > 3).any() else 'control', axis=1)# july-changes
filtered_stuff=annot_5_tissues.groupby('gene').filter(lambda group: 'outlier' in set(group.Y)).reset_index(drop=True)
no_bnd_inv_filtered=filtered_stuff[(~filtered_stuff.SV.str.contains('BND'))&(~filtered_stuff.SV.str.contains('INV'))]
no_bnd_inv_filtered=no_bnd_inv_filtered.loc[:, (no_bnd_inv_filtered != no_bnd_inv_filtered.iloc[0]).any()] 
SV_count=var_matrix.iloc[:,:-1].gt(0).sum(axis=1).reset_index().rename({'Unnamed: 0':'SV',0:'SV_ac'},axis=1)

no_bnd_inv_filtered = no_bnd_inv_filtered.merge(SV_count,how='left',on='SV')
no_bnd_inv_filtered['CN_change'] = no_bnd_inv_filtered.apply(lambda x: cn_change(x['type']),axis=1)

dummified = pd.get_dummies(no_bnd_inv_filtered,columns = ['type','CN_change']).rename({'Ind':'SubjectID','gene':'GeneName'},axis=1).drop_duplicates()

args = {}
with open('starting_annotations/tissue_gene_level_methods.tsv','r') as handle:
    for line in handle:
        feature,method = line.rstrip().split('\t')
        args[feature] = method

args['SV_ac'] = 'min'
args = {key:value for key,value in args.items() if key in dummified.columns}

gene_level_annotations = dummified.groupby(['SubjectID','GeneName']).agg(args).reset_index()
gene_level_sv_count = dummified.groupby(['SubjectID','GeneName']).agg(num_rare_var=('SV','nunique')).reset_index()
gene_level_annotations=gene_level_annotations.merge(gene_level_sv_count,how='left',on=['SubjectID','GeneName'])

#gene_level_annotations['TE_pvalues']
tissue_columns = [i for i in gene_level_annotations.columns if 'Z:' in i]
tissue_z_to_p_names = {i:f'P:{i.split(":")[1]}' for i in tissue_columns}
gene_level_annotations[tissue_columns]=gene_level_annotations[tissue_columns].apply(Zscores_to_Pvalues,axis=0)

gene_level_annotations.rename(columns=tissue_z_to_p_names,inplace=True)

pre_n2pair=gene_level_annotations.drop(['Y'],axis=1)
grouping=dummified.groupby(['SubjectID','GeneName']).agg(SVs=('SV',lambda x: ','.join(np.sort(x.unique()))),
                                                         Y=('Y','first'),).reset_index()
grouping=grouping.merge(gene_indiv_set,how="right",on=["SubjectID","GeneName"])
out=[]
n2_id=0
for index,df in grouping.groupby(['SVs','GeneName']):
    if df.shape[0]>1:
        n2pair=df.sample(n=2)
        n2pair['N2pair']=n2_id
        n2_id+=1
        out.append(n2pair)
    else:
        df['N2pair']='NA'
        out.append(df)
        
training_testing=pd.concat(out)

train_test_split=training_testing.drop(['SVs','Y'],axis=1)

test_debug=gene_level_annotations.merge(train_test_split,how='right',on=['SubjectID','GeneName'])

test_debug['type_INS'] = test_debug[['type_ALU', 'type_LINE1', 'type_SVA']].\
apply(lambda x: 1 if sum(x.to_numpy()) > 0 else 0,axis=1)
test_debug['type_DEL'] = test_debug[['type_MEI', 'type_DEL']].\
apply(lambda x: 1 if sum(x.to_numpy()) > 0 else 0,axis=1)
test_debug=test_debug.drop(['type_ALU', 'type_LINE1', 'type_SVA','type_MEI'],axis=1)


test_debug.to_csv('tissue_Watershed_annotation_july21.csv',sep='\t',header=True,index=False)
