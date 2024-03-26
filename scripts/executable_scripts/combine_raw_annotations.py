import pandas as pd
import numpy as np
import scipy.stats as stats
import argparse
from sv_utils import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='filter and combine to data set.')
    #parser.add_argument('--sv-info',type=str,required=True,metavar='[pipeline input bed]',
    #                help='the bed file containing additional info: SV,SVTYPE')
    parser.add_argument('--maf-frame',type=str,required=True,metavar='[maf tsv]',
                    help='maf information')
    parser.add_argument('--genotype',type=str,metavar='[pipeline input allele file]',required=True,
                    help='allele file see formats')
    parser.add_argument('--gene-sv',type=str,metavar='[gene sv file]',required=True,
                    help='path to output the genotypes to. optional if running without genotypes.')
    parser.add_argument('--outlier-file',type=str,required=True,metavar='[expression file]',
                    help='extract only samples with genotype filter in filters(ie. "PASS . MATCH_1KGP") provided')
    parser.add_argument('--additional-annotations',type=str,nargs='+',metavar='[annot1 annot2 annot3...]',
                    default=[],help='a space separated list of annotation files to incorporate.')
    parser.add_argument('--out-raw-annotations',type=str,required=True,metavar='[annotation file with no additional annotations]',
                    help='path to annotation file, later will append additional annotation columns to this one.')
    

    args = parser.parse_args()
    # input block.
    sv_file=args.sv_info
    maf_file=args.maf_frame
    genotype_file=args.genotype
    gene_sv_file=args.gene_sv
    outlier_file=args.outlier_file
    annotation_file_list=args.additional_annotations
    output=args.out_raw_annotations
    #svs = pd.read_csv(sv_file,sep='\t',header=None)
    #mafs = pd.read_csv(maf_file,sep='\t')
    genotypes = pd.read_csv(genotype_file,sep='\t')
    gene_sv = pd.read_csv(gene_sv_file,sep='\t',header=None,names=['chr','start','end','SVid','SVTYPE','Gene'])
    medianZ_outliers = pd.read_csv(outlier_file,sep='\t')

    # check if variants are flipped:
    # filter 1, by qual filter
    filtered_var = genotypes[genotypes.Allele>=0]

    # filter 2, by not using the filter function directly, we avoided the loop. 
    rare_status = filtered_var.groupby('SV').apply(is_rare)
    # then calculate the condition first 
    # finally use loc/iloc to get cols and rows of interest for maximal performance
    rare_index = rare_status[rare_status==True].index
    rare_SV_cond = filtered_var.SV.isin(rare_index)
    rare_var = filtered_var.loc[rare_SV_cond]

    # filter 3, flip variants if MAF flipped...
    flipped_status = rare_var.groupby('SV').apply(is_flipped)
    flipped_index = flipped_status[flipped_status==True].index
    flipped_cond = rare_var.SV.isin(flipped_index)
    rare_var.SVTYPE = collapse_types(rare_var.SVTYPE)
    rare_var.loc[flipped_cond,'Allele'] = 2-rare_var.loc[flipped_cond,'Allele']
    rare_var.loc[flipped_cond,'SVTYPE'] = flipped_type(rare_var.loc[flipped_cond,'SVTYPE'])
    rare_var['SVid'] = [i.split('.')[0] for i in rare_var.SV]


    # compute allele frequency here
    maf_frame = rare_var.groupby('SV').apply(maf).reset_index(name='af')
    # filter again for rare cnv
    maf_frame['SVid'] = maf_frame['SV'].str.split('.').str[0]
    maf_frame_collapse_maf=maf_frame.merge(maf_frame.groupby('SVid')['af'].sum().reset_index(),on='SVid',suffixes=('','_collapsed'))
    true_rare_maf_frame = maf_frame_collapse_maf[maf_frame_collapse_maf.af_collapsed<0.01]

    # next step is to filter to nonzero and high qual alleles.
    true_rare_var = rare_var[rare_var.SV.isin(true_rare_maf_frame.SV)]
    true_rare_alleles = true_rare_var[(true_rare_var['Allele']>0)].copy()

    # then filter to those that are near by genes,
    SV_Ind_Gene_noExpn = true_rare_alleles.merge(gene_sv[['Gene','SVid']],on=['SVid'],how='left')

    # SV->Ind->Gene(further see which one has measured expression)
    SV_Ind_Gene = SV_Ind_Gene_noExpn.merge(medianZ_outliers.rename(columns={'Ind':'SUBJID'}),on=['SUBJID','Gene'],how='inner')


    # filter to no breakend
    cond_noBNDINV=(SV_Ind_Gene['SVTYPE'] !='BND')
    noBNDINV_SV_Ind_Gene = SV_Ind_Gene[cond_noBNDINV]
    # filter to make sure at least 1 individual is outlier. 
    gene_keep_condition_by_Y = noBNDINV_SV_Ind_Gene.groupby('Gene').apply(lambda x: True if np.count_nonzero(x['Y'].to_numpy()=='outlier') >0 else False)
    gene_kept = gene_keep_condition_by_Y[gene_keep_condition_by_Y==True].index.unique()
    index_kept = noBNDINV_SV_Ind_Gene.Gene.isin(gene_kept)
    filtered_total_data = noBNDINV_SV_Ind_Gene.loc[index_kept,:]
    #filtered_total_data.to_csv(output,sep='\t',header=True,index=False)

    # combine with other annotations. make sure they have SV gene field. 
    annots = []
    for annot_file in annotation_file_list:
        annots.append(pd.read_csv(annot_file,sep='\t'))
    merged_annotations=functools.reduce(lambda left,right: pd.merge(left,right,how='left',on=['Gene','SV']),annots)
    uncollapsed=merged_annotations.rename(columns={"Gene":"GeneName","Ind":"SubjectID"})
    uncollapsed.to_csv(output,sep='\t',header=True,index=False)


