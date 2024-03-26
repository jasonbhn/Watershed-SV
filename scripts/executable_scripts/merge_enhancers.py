import pandas as pd
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='merge enhancer annotations')
    parser.add_argument('--enhancers',type=str,required=True,metavar='[enhancer file]',
                    help='enhancer file')
    parser.add_argument('--primary-cell-list',type=str,required=True,metavar='[list of names of primary cell types]',
                    help='list of names of primary cell types. ')
    parser.add_argument('--out-merged-enhancers',type=str,required=True,metavar='[merged enhancers primary cell]',
                    help='output merged enhancer cell types')
    args = parser.parse_args()
    # input argument variables
    enhancers = args.enhancers
    primary_cell_list = args.primary_cell_list
    out_merged_enhancers = args.out_merged_enhancers


    multitissue_enhancers=pd.read_csv(enhancers)
    primary_cells=[]
    with open(primary_cell_list,'r')as handle:
        for line in handle:
            primary_cells.append(line.rstrip())
    primary_enhancers=multitissue_enhancers[['Enh_peaks']+primary_cells]

    primary_enhancers=primary_enhancers.set_index('Enh_peaks')
    primary_enhancers=primary_enhancers[(primary_enhancers.T!=0).any()]
    primary_enhancers['num_cell_types']=primary_enhancers.apply(np.count_nonzero,axis=1)
    primary_enhancers['chrom']=primary_enhancers.index.str.split(':').str[0]
    primary_enhancers['start_end_str']=primary_enhancers.index.str.split(':').str[1]

    primary_enhancers['start']=primary_enhancers.start_end_str.str.split('-').str[0].astype('int')
    primary_enhancers['end']=primary_enhancers.start_end_str.str.split('-').str[1].astype('int')
    primary_enhancers[['chrom','start','end','num_cell_types']].reset_index().drop('Enh_peaks',axis=1).to_csv(out_merged_enhancers,sep='\t',header=False,index=False)



