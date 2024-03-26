import pandas as pd
import numpy as np
import pyranges as pr
import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='tools to clean up, impute watershed input')
    parser.add_argument('--cluster',type=str,required=True,metavar='[cluster file]',
                    help='parse this cohort cluster file')
    parser.add_argument('--work-dir',type=str,required=True,metavar='[where to output filtered junction files]')
    parser.add_argument('--filtered-junclist',type=str,required=True,metavar='[where to output fitlered junc file lists]',
                    help='set of intron bed coord to keep.')
    parser.add_argument('--junclist',type=str,required=True,metavar='[list containing the unfiltered junctions]',
                    help='unfiltered intron clusters')
    args = parser.parse_args()
unfiltered_clusters = args.cluster
workdir = args.work_dir
juncListFile = args.junclist
filtered_JuncList = args.filtered_junclist
splice_scan = pd.read_csv(unfiltered_clusters,
                          sep=' ')
sparse_cluster_index = np.where(np.any(splice_scan.to_numpy() >= 15,axis=1)==True)[0]
reindexed_splice_scan = splice_scan.reset_index()
reindexed_splice_scan_removeSparseJunc = reindexed_splice_scan.loc[sparse_cluster_index]
reindexed_splice_scan_removeSparseJunc[['chrom','start','end','cluster']] = reindexed_splice_scan_removeSparseJunc['index'].str.split(':',expand=True)
reindexed_splice_scan_removeSparseJunc[['start','end']] = reindexed_splice_scan_removeSparseJunc[['start','end']].astype(int)
# filter with pyranges with dask parallel process. 
with open (juncListFile, 'r') as handle:
    juncFiles = [juncFile.rstrip() for juncFile in handle.readlines()]
filteredJuncList = []
for juncFile in juncFiles:
    filteredJuncFile = f'{os.path.join(workdir,os.path.basename(juncFile,'junc'))}filtered.junc'
    data = pd.read_csv(juncFile,
    sep='\t',
    header=None,
    names=['chrom','chromStart','chromEnd','name','score',
    'strand','thickStart','thickEnd','itemRgb','blockCount',
    'blockSizes','blockStarts'])
    data[['StartSize','EndSize']] = data['blockSizes'].str.split(',',expand=True).astype(int)
    data['start'] = data['StartSize'] + data['chromStart']
    data['end'] = data['chromEnd'] - data['EndSize'] + 1
    filtered_data = data.\
        merge(reindexed_splice_scan_removeSparseJunc[['chrom','start','end']],
        on=['chrom','start','end'],how='inner')
    filtered_data[['chrom','chromStart','chromEnd','name','score',
    'strand','thickStart','thickEnd','itemRgb','blockCount',
    'blockSizes','blockStarts']].to_csv(filteredJuncFile,sep='\t',header=False,index=False)
    filteredJuncList.append(filteredJuncFile)
with open(filtered_JuncList,'w') as handle:
    for i in filteredJuncList:
        handle.write(f'{i}\n')
        
