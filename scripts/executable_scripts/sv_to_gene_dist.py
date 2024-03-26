import pandas as pd
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gene to SV distance')
    parser.add_argument('--gene-sv',type=str,required=True,metavar='[gene SV pairing]',
                    help='gene sv pair file ')
    parser.add_argument('--gene',type=str,required=True,metavar='[gene info]',
                    help='gene file')
    parser.add_argument('--flank',type=int,required=True,metavar='[flank used for SV annotations]',
                    help='flank distance')
    parser.add_argument('--in-gene-tss',type=str,required=True,metavar='[input for SV-tss pair annotations]',
                    help='path to input file to be created for SV-tss annotations')
    parser.add_argument('--in-gene-tes',type=str,required=True,metavar='[input for SV-tes pair annotations]',
                    help='path to input file to be created for SV-tes annotations')
    parser.add_argument('--out-gene-sv-dist',type=str,required=True,metavar='distance of sv to gene',
                    help='output the gene sv distance info. ')
    args = parser.parse_args()
    # input argument variables
    gene_sv = args.gene_sv
    gene = args.gene
    flank = args.flank
    in_gene_tss = args.in_gene_tss
    in_gene_tes = args.in_gene_tes
    out_gene_sv_dist = args.out_gene_sv_dist
    def determine_noncoding_distance_to_tss(tss,olstart,olend,strand,mode='tss'):
        if mode=='tss':
            if strand == '+':
                if olend > tss:
                    return 0
                else:
                    return tss - olend
            if strand == '-':
                if olstart < tss: 
                    return 0
                else:
                    return olstart - tss
        else:
            if strand == '+':
                if olstart < tss:
                    return 0
                else:
                    return olstart - tss
            if strand == '-':
                if olend > tss: 
                    return 0
                else:
                    return tss - olend
    def determine_noncoding_percentage(tss,olstart,olend,strand,flanking_distance,mode='tss'):
        if mode=='tss':
            if strand == '+':
                if olend > tss:
                    return 0
                else:
                    return olend - max(tss - flanking_distance,olstart)
            if strand == '-':
                if olstart < tss: 
                    return 0
                else:
                    return min(tss + flanking_distance,olend) - olstart
        else:
            if strand == '+':
                if olstart < tss:
                    return 0
                else:
                    return min(tss + flanking_distance,olend) - olstart
            if strand == '-':
                if olend > tss: 
                    return 0
                else:
                    return olend - max(tss - flanking_distance,olstart)
    gene_sv=pd.read_csv(gene_sv,sep='\t',header=None,names=['chrom','olstart','olend','SV','SVTYPE','gene'])
    gene_strand=pd.read_csv(gene,sep='\t',header=None,names=['chrom','gene_start','gene_end','gene','biotype','strand'])
    gene_tss=pd.read_csv(in_gene_tss,sep='\t',header=None,names=['_1','tss','_2','gene']).drop(['_1','_2'],axis=1)
    gene_tes=pd.read_csv(in_gene_tes,sep='\t',header=None,names=['_1','tes','_2','gene']).drop(['_1','_2'],axis=1)
    genes_tss_tes=gene_sv.merge(gene_tss,how='left',on='gene').\
        merge(gene_tes,how='left',on='gene').\
            merge(gene_strand[['gene','strand']],how='left',on='gene')
    genes_tss_tes['nc_dist2tss'] = genes_tss_tes.apply(lambda x: determine_noncoding_distance_to_tss(x.tss,x.olstart,x.olend,x.strand),axis=1)
    genes_tss_tes['tss_flank_impacted'] = genes_tss_tes.apply(lambda x: determine_noncoding_percentage(x.tss,x.olstart,x.olend,x.strand,flank),axis=1)
    #genes_tss_tes['dist2tss'] = genes_tss_tes.apply(lambda x: 0 if (x.olstart <x.tss and x.olend >x.tss) else min([abs(x.olstart-x.tss),abs(x.olend-x.tss)]),axis=1)
    genes_tss_tes['nc_dist2tes'] = genes_tss_tes.apply(lambda x: determine_noncoding_distance_to_tss(x.tes,x.olstart,x.olend,x.strand,mode='-'),axis=1)
    genes_tss_tes['tes_flank_impacted'] = genes_tss_tes.apply(lambda x: determine_noncoding_percentage(x.tes,x.olstart,x.olend,x.strand,flank,mode='-'),axis=1)
    #genes_tss_tes['dist2tes'] = genes_tss_tes.apply(lambda x: 0 if (x.olstart <x.tes and x.olend >x.tes) else min([abs(x.olstart-x.tes),abs(x.olend-x.tes)]),axis=1)
    genes_tss_tes['nc_dist2gene'] = genes_tss_tes['nc_dist2tss'] + genes_tss_tes['nc_dist2tes']
    genes_tss_tes['upstream_noncoding_variant'] = np.where(genes_tss_tes['nc_dist2tss']>0,1,0)
    genes_tss_tes['downstream_noncoding_variant'] = np.where(genes_tss_tes['nc_dist2tes']>0,1,0)
    genes_tss_tes[['gene','SV','nc_dist2gene','upstream_noncoding_variant','downstream_noncoding_variant','tss_flank_impacted','tes_flank_impacted']].rename(columns={'gene':'Gene'}).to_csv(out_gene_sv_dist,sep='\t',header=True,index=False)
