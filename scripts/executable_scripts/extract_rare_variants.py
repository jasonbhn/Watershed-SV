#!/usr/bin/env python3

import sys
from pysam import VariantFile
import pandas as pd
import numpy as np
from sv_utils import *
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='VCF filtering based on allele frequency field')
    parser.add_argument('--vcf',type=str,required=True,metavar='[input rare SV VCF file]',
                    help='select [rare_SV].vcf to parse')
    parser.add_argument('--lifted-coord',type=str,metavar='[liftover coordinates]',
                    default="None",
                    help='provide liftovered bed coordinate of SV format is: [CHR|START|END|SVID], if vcf in 38, ignore this argument')
    parser.add_argument('--metadata',type=str,metavar='[metadata for sample filtering]',
                    default="None",
                    help='provide metadata for filtering samples in population vcf')
    parser.add_argument('--out-genotype',type=str,metavar='[output genotypes]',
                    default="None",
                    help='path to output the genotypes to. optional if running without genotypes.')
    parser.add_argument('--genotype-filters',type=str,nargs="+",required=True,metavar='[genotype filters]',
                    default=["PASS",".","MATCH_1KGP"],
                    help='extract only samples with genotype filter in filters(ie. "PASS . MATCH_1KGP") provided')
    parser.add_argument('--out-annotsv',type=str,required=True,metavar='[input for downstream annotSV]',
                    help='path to input file to be created for annotSV')
    parser.add_argument('--out-generic',type=str,required=True,metavar='[annotation file with no additional annotations]',
                    help='path to annotation file, later will append additional annotation columns to this one.')
    parser.add_argument('--out-maf',type=str,required=False,metavar='[maf file]',
                    help='path to maf file, later will append additional annotation columns to this one.')
    parser.add_argument('--extract-genotype',action=argparse.BooleanOptionalAction,metavar='[--extract-genotype|--no-extract-genotype]',
                    help='extract genotypes or not.')
    parser.add_argument('--filter-ethnicity',action=argparse.BooleanOptionalAction,metavar='[--filter-ethnicity|--no-filter-ethnicity]',
                    help='filter ethnicity or not.')
    parser.add_argument('--infer-rareness',action=argparse.BooleanOptionalAction,metavar='[--extract-genotype|--no-extract-genotype]',
                    help='whether to infer rare vars from vcf or not.')
    

    args = parser.parse_args()
    # input argument variables
    vcf = args.vcf # VCF file of interest. If to use with Watershed, must be rare checked on your own!
    lifted_coord = args.lifted_coord # if no --lifted-coord, then is string "None"
    out_genotype = args.out_genotype # if no --out-genotype, then is string "None"
    genotype_filters = args.genotype_filters # this is a list of filters
    out_annotsv = args.out_annotsv
    output = args.out_generic
    out_maf = args.out_maf
    extract_genotype = args.extract_genotype
    infer_rareness = args.infer_rareness
    filter_ethnicity = args.filter_ethnicity
    metadata_path = args.metadata
    # extract genotypes into a file. 
    if lifted_coord == "None":
        print('hg38 files are provided, if not, please provide lifted bed file in hg38')
        personal_genotypes_list,bed_df = extract_genotype_tuples(vcf,\
        filters_pass=genotype_filters,extract_genotypes=extract_genotype,\
        filter_ethnicity=filter_ethnicity,metadata_path=metadata_path)
    else:
        print('changing coordinates to lifted coord in hg38')        
        personal_genotypes_list,unlifted_bed_df = extract_genotype_tuples(vcf,\
        filters_pass=genotype_filters,extract_genotypes=extract_genotype,\
        filter_ethnicity=filter_ethnicity,metadata_path=metadata_path)
        hg38_coord = pd.read_csv(lifted_coord,sep='\t',header=None,names=['chrom','start','end','SV','MapQ'])
        if hg38_coord['SV'].str.contains('CNV').any():
            cnv_coord=hg38_coord[hg38_coord['SV'].str.contains('CNV')]
            dup_cnv_coord = cnv_coord.copy()
            dup_cnv_coord['SV'] = dup_cnv_coord['SV'] + '.DUP_CNV'
            del_cnv_coord = cnv_coord.copy()
            del_cnv_coord['SV'] = del_cnv_coord['SV'] + '.DEL_CNV'
            complete_hg38_coord = pd.concat([hg38_coord,dup_cnv_coord,del_cnv_coord])
        info_cols = [i for i in unlifted_bed_df.columns if i not in ['chrom','start','end']]
        bed_df = unlifted_bed_df[info_cols].merge(complete_hg38_coord,on='SV',how='inner')[['chrom','start','end']+info_cols]
        # done lifting the coord. 
    # now start generate the three important outputs. 
    var_bed = bed_df[['chrom','start','end','SV','SVTYPE']].copy()
    var_bed['collaped_SVTYPE'] = collapse_types(var_bed['SVTYPE'].tolist())
    var_bed['VEP_SVTYPE'] = collapse_types(var_bed['SVTYPE'].tolist())
    var_bed['annotSV_SVTYPE'] = svtype_annotsv(var_bed['SVTYPE'].tolist())
    var_bed['end'] = np.where(var_bed.start==var_bed.end,var_bed.start+1,var_bed.end)
    # when extract_genotype == True, output genotypes:
    if extract_genotype:
        genotypes = pd.concat([pd.DataFrame(data=personal_genotypes) for personal_genotypes in personal_genotypes_list])
        genotypes.columns = ['SUBJID','SV','SVTYPE','Allele','CN']
    
    if infer_rareness:
        # check if variants are flipped:
        # filter 1, by qual filter
        genotypes['Allele']= genotypes['Allele'].astype('int')
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
        #true_rare_maf_frame = maf_frame_collapse_maf[maf_frame_collapse_maf.af_collapsed<0.01]
        # write true rare var's maf to a dataframe. 
        maf_frame.to_csv(out_maf,sep='\t',header=True,index=False)
        # next step is to filter to nonzero and high qual alleles.
        true_rare_var = rare_var[rare_var.SV.isin(maf_frame.SV)]
        true_rare_alleles = true_rare_var[(true_rare_var['Allele']>0)].copy()
        allvar_bed = var_bed[(var_bed.start>=0)&(var_bed.end>0)][['chrom','start','end','SV','collaped_SVTYPE']]
        allvar_bed[allvar_bed['SV'].isin(true_rare_alleles[true_rare_alleles.SVTYPE!='BND'].SV.unique())].to_csv(output,index=False,header=False,sep='\t')
        allvar_annotsv_bed = var_bed[(var_bed.start>=0)&(var_bed.end>0)][['chrom','start','end','SV','annotSV_SVTYPE']]
        allvar_annotsv_bed[allvar_annotsv_bed['SV'].isin(true_rare_alleles[true_rare_alleles.SVTYPE!='BND'].SVid.unique())].to_csv(out_annotsv,index=False,header=False,sep='\t')
        true_rare_alleles = true_rare_var[(true_rare_var['Allele']>0)].copy()
        true_rare_alleles.to_csv(out_genotype,sep='\t',header=True,index=False)

    else: 
        genotypes.to_csv(out_genotype,sep='\t',header=True,index=False)
        var_bed[(var_bed.start>=0)&(var_bed.end>0)][['chrom','start','end','SV','collaped_SVTYPE']].to_csv(output,index=False,header=False,sep='\t')
        var_bed[(var_bed.start>=0)&(var_bed.end>0)][['chrom','start','end','SV','annotSV_SVTYPE']].to_csv(out_annotsv,index=False,header=False,sep='\t')
