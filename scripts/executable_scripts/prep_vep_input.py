#!/usr/bin/env python3

import pandas as pd
import sys
def svtype_vep(col):
    if col == 'INS':
        return 'INS'
    elif col == 'DEL' or col == 'DEL_CNV':
        return 'DEL'
    elif col == 'DUP':
        return 'DUP'
    elif col == 'DUP_CNV':
        return 'TDUP'
    elif col =='INV':
        return 'DEL'
gene_sv=pd.read_csv(sys.argv[1],
                    sep='\t',
                    header=None,
                    names=['chrom','start','end','SV','SVTYPE','Gene'])
# get vep input
vep_upload=gene_sv[['chrom','start','end','SV','SVTYPE','Gene']]
vep_upload['id'] = vep_upload.SV + ":" + vep_upload.Gene
vep_upload['type'] = vep_upload.apply(lambda x: svtype_vep(x['SVTYPE']),axis=1)
vep_upload['strand'] = '+'
vep_upload = vep_upload.drop(['SV','Gene','SVTYPE'],axis=1)
vep_upload[['chrom','start','end','type','strand','id']].to_csv(sys.argv[2],header=False,index=False,sep='\t')
