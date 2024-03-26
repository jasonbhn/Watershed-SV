import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels as sm
import seaborn as sns
import pickle
from pandarallel import pandarallel
from matplotlib import pyplot as plt
from sklearn.impute import KNNImputer
from statsmodels.stats.multitest import fdrcorrection
pandarallel.initialize(nb_workers=24,progress_bar=True)
def impute(df,k=200):
    imputer = KNNImputer(n_neighbors=k)
    data = df.iloc[:,2:].to_numpy().T
    id_cols = df.iloc[:,1:2]
    filled_data = pd.DataFrame(data=imputer.fit_transform(data).T,columns=df.iloc[:,2:].columns,index=df.iloc[:,2:].index)
    return pd.concat([id_cols, filled_data],axis=1).set_index('tissue')

multitissue_gene_indiv_zscore = pd.read_csv('test_correlation_joint_data_v3.csv',sep='\t')
# filter tissue first, ind second. 
filter_by_tissue = multitissue_gene_indiv_zscore[(multitissue_gene_indiv_zscore.count(axis=1)-2)/(multitissue_gene_indiv_zscore.shape[1]-2)>0.25]
filter_by_tissue_ind = filter_by_tissue.loc[:,(filter_by_tissue.count()/filter_by_tissue.shape[0])>=0.25]

imputed_data = filter_by_tissue_ind.groupby('Gene').parallel_apply(lambda x: impute(x)).reset_index()
imputed_data.to_csv('imputed_correlational_outlier_data_pandarallel.csv',index=False,sep='\t')
