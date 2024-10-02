#!/usr/bin/env python3

import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gene tp SV distance')
    parser.add_argument('--gene-sv-cpg',type=str,required=True,metavar='[gene SV cpg pairing]',
                    help='gene sv cpg file')
    parser.add_argument('--out-gene-sv-cpg',type=str,required=True,metavar='collapsed sv to gene cpg',
                    help='output the gene sv cpg info. ')
    args = parser.parse_args()
    # input argument variables
    sv_gene_cpg = args.gene_sv_cpg
    out_gene_sv_cpg = args.out_gene_sv_cpg
    sv_gene_cpg = pd.read_csv(sv_gene_cpg,sep='\t',header=None)
    sv_gene_cpg.columns=['chrom','start','end','SV','SVTYPE','Gene','cpgchrom','cpgstart','cpgend','cpgpct']
    maxcpg=sv_gene_cpg.groupby(['Gene','SV']).agg(max_cpgpct=('cpgpct','max')).reset_index()
    maxcpg.to_csv(out_gene_sv_cpg,sep='\t',header=True,index=None)
