from scipy import stats
import pandas as pd
import numpy as np 
import glob
import argparse
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='gene tp SV distance')
    parser.add_argument('--tissue-expn-dir',type=str,required=True,metavar='[directory containing tissue expression]',
                    help='gene sv cpg file')
    parser.add_argument('--mode',default='median',choices=['tissue-specific', 'median'],metavar='[whether to keep tissue specific outliers, or collapse to median outlier]',
                    help='which outlier to call')
    parser.add_arl
    args = parser.parse_args()
    input_files=args.tissue_expn_dir

    i=0
    tissue_list=[]
    for input_file in input_files:
        tissue=input_file.split('/')[-1].split('.')[0]
        ste_unstack=pd.read_csv(input_file,sep='\t')
        ste = ste_unstack.set_index('Id').stack().reset_index().rename(columns={'Id':'Gene','level_1':'Ind',0:'Zscore'})
        ste=ste.rename(columns={'Zscore':f'Z:{tissue}'})
        tissue_list.append(ste)
        i+=1
        print(f'added with {i} tissues')
    primer_frame=tissue_list[0]
    for i in range(1,len(tissue_list)):
        primer_frame=primer_frame.merge(tissue_list[i],on=['Gene','Ind'],how='outer')
        print(f'merged with {i} tissues')
    #write it to file first
    #primer_frame.to_csv(snakemake.output['multiZ'],sep='\t',header=True,index=False)
    indexed_primer_frame=primer_frame.set_index(['Gene','Ind'])
    num_tissue_filter=indexed_primer_frame.count(axis=1)>=5
    MedZ_frame = indexed_primer_frame.loc[num_tissue_filter,:].median(axis=1).reset_index().rename(columns={0:'MedZ'})
    MedZ_frame['Y'] = np.where(MedZ_frame['MedZ'].abs() > snakemake.params['medz_cutoff'], 'outlier','control')
    # now groupby SubjectID to eliminate global outliers.
    oc_prop=MedZ_frame.groupby('Ind')['Y'].value_counts(normalize=True)
    q3=oc_prop[:, 'outlier'].quantile(0.75)
    q1=oc_prop[:, 'outlier'].quantile(0.25)
    global_outlier_threshold=1.5*(q3-q1) + q3
    drop_SubjectIDs = oc_prop[:, 'outlier'].index[oc_prop[:, 'outlier'] > global_outlier_threshold]
    keep_SubjectIDs = (~MedZ_frame['Ind'].isin(drop_SubjectIDs))
    noGlobalOutlier_MedZ_frame = MedZ_frame[keep_SubjectIDs]
    noGlobalOutlier_MedZ_frame.to_csv(snakemake.output['medZ'],header=True,index=False,sep='\t')

