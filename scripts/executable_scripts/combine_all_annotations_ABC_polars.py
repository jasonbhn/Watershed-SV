#!/usr/bin/env python3

from pysam import VariantFile
from sv_utils import *
import pandas as pd
import polars as pl
import polars.selectors as cs
import functools
import os
import numpy as np
import argparse

def extract_MAF(vcf_path,field):
    additional_SV_info = []
    vcf_read = VariantFile(vcf_path,'r')
    for rec in vcf_read.fetch():
        additional_SV_info.append([rec.id,rec.info[field]])
    return pd.DataFrame(additional_SV_info,columns=['SV',field])
def extract_length(vcf_path,field='SVLEN'):
    additional_SV_info = []
    vcf_read = VariantFile(vcf_path,'r')
    for rec in vcf_read.fetch():
        if rec.info['SVTYPE']!='BND':
            if rec.info['SVTYPE']=='LINE1':
                if isinstance(rec.info['SVLEN'],int):
                    if rec.info['SVLEN']!= -1:
                        additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN']))])
                    else:
                        additional_SV_info.append([rec.id,abs(int(rec.info['MEINFO'][2]))])
                else: 
                    if rec.info['SVLEN'][0]!= -1:
                        additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN'][0]))])
                    else:
                        additional_SV_info.append([rec.id,abs(int(rec.info['MEINFO'][2]))])
            elif rec.info['SVTYPE']== 'MEI':
                if isinstance(rec.info['SVLEN'],int):
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN']))])
                else:
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN'][0]))])
            elif rec.info['SVTYPE']== 'SVA':
                if isinstance(rec.info['SVLEN'],int):
                    svlen = abs(int(rec.info['SVLEN']))
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN']))])
                else:
                    svlen = abs(int(rec.info['SVLEN'][0]))
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN'][0]))])
            elif rec.info['SVTYPE']== 'ALU':
                if isinstance(rec.info['SVLEN'],int):
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN']))])
                else:
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN'][0]))])
            elif rec.info['SVTYPE']=='BND':
                additional_SV_info.append([rec.id,-1])
            elif 'CNV' in rec.info['SVTYPE']:
                if 'SVLEN' in rec.info:
                    if isinstance(rec.info['SVLEN'],int):
                        additional_SV_info.append([rec.id+'.DUP_CNV',abs(int(rec.info['SVLEN']))])
                        additional_SV_info.append([rec.id+'.DEL_CNV',abs(int(rec.info['SVLEN']))])
                    else:
                        additional_SV_info.append([rec.id+'.DUP_CNV',abs(int(rec.info['SVLEN'][0]))])
                        additional_SV_info.append([rec.id+'.DEL_CNV',abs(int(rec.info['SVLEN'][0]))])
                else:
                    additional_SV_info.append([rec.id+'.DUP_CNV',abs(rec.stop-rec.pos)])
                    additional_SV_info.append([rec.id+'.DEL_CNV',abs(rec.stop-rec.pos)])
            elif rec.info['SVTYPE']=='DEL':
                if isinstance(rec.info['SVLEN'],int):
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN']))])
                else:
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN'][0]))])
            else:
                if isinstance(rec.info['SVLEN'],int):
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN']))])
                else:
                    additional_SV_info.append([rec.id,abs(int(rec.info['SVLEN'][0]))])

    return pd.DataFrame(additional_SV_info,columns=['SV',field])



# several different modes 
# 1. collapse one single annotation set. 
#   a. maf needs to be calculated here:
#       i. user supply custom MAF, required format is SVid -> MAF
#       ii. user point out which info field name, typically 'MAF', program extracts from vcf
#   b. length needs to be calculated here:
#       i. user supply custom length, required format is SVid -> length
#       ii. user point to info field name, typically 'SVLEN', program extracts from vcf
#   c. CN calculation is different based on dataset (ie, long /short/ cnv)
#       i. user supply custom CN, required format is SVid -> CN (note this can be tandem expansion)
#       ii. user point to info field name, typically 'CN', program extracts from vcf
#       !!!! note that CN are SubjectID->SV specific, so the file needs to have 3 columns. 
#       !!!! note the proragm will impute CN if -1 regardless of i/ii. 
#   d. collapse to gene level annotation or keep var: a flag 
#
# 2. combine two set to trim down the features based on presence abscence
#   a. mode i: generate eval_watershed input, split like before
#   b. mode ii: combine training + testing, trim and collapse, generate two outfiles
#   c. output files prefix, mode i, 1 file, mode ii 2 files. 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='tools to clean up, impute watershed input')
    parser.add_argument('--vcf',type=str,required=True,metavar='[input rare SV VCF file]',
                    help='select [rare_SV].vcf to parse')
    parser.add_argument('--genotypes',type=str,required=True,metavar='[input genotype VCF file]',
                    help='select [rare_SV].genotype to parse')
    parser.add_argument('--genes',type=str,required=True,metavar='[input gene list file]',
                    help='select genes to parse')
    parser.add_argument('--gene-sv',type=str,required=True,metavar='[input sv-gene pair file]',
                    help='select sv gene pairs to parse')
    parser.add_argument('--annotation-dir',type=str,required=True,metavar='[where the annotations are generater]',
                    help='directory name of annotations, usually something.../intermediates, use anything before intermediates')
    parser.add_argument('--outfile',type=str,required=True,metavar='[where to output annotation]',
                    help='where to output single set annotation. ')
    parser.add_argument('--expressions',type=str,required=True,metavar='[input gene expression file]',
                    help='select gene expression to parse')
    parser.add_argument('--expression-field',type=str,required=True,metavar='[column name of zscores]',
                    help='select gene expression column to parse, it could be the PREFIX of columns, like TE_pvalues_Blood, TE_pvalues_Muscle...')
    parser.add_argument('--expression-id-field',type=str,required=True,metavar='[column name of each sample]',
                    help='select sample column to parse')
    parser.add_argument('--maf-mode',required=True,choices=['upload', 'extract'],metavar='[choose how to get MAF]',
                    help='either upload custom MAF SV->MAF tsv, or extract from info from vcf')
    parser.add_argument('--maf-file',default='None',type=str,metavar='[if --maf-mode upload, give file name]',
                    help='name of tsv with SVid from vcf, and a maf column')
    parser.add_argument('--maf-field',default='Max_AF',type=str,metavar='[if --maf-mode extract, give field name]',
                    help='name of info field from vcf. default follows SVAFotate.')
    parser.add_argument('--length-mode',required=True,choices=['upload-SV', 'upload-VNTR', 'extract'],metavar='[choose how to get length]',
                    help='either upload custom MAF SV->length tsv, or extract from info from vcf')
    parser.add_argument('--length-file',default='None',type=str,metavar='[if --length-mode upload, give file name]',
                    help='name of tsv with SVid from vcf, and a length column')
    parser.add_argument('--length-field',default='SVLEN',type=str,metavar='[if --length-mode extract, give field name]',
                    help='name of info field from vcf. default follows vcf 4.2.')
    parser.add_argument('--CN-mode',required=True,choices=['upload', 'extract'],metavar='[choose how to get CN]',
                    help='either upload custom MAF SV->CN tsv, or extract from info from vcf')
    parser.add_argument('--CN-file',default='None',type=str,metavar='[if --CN-mode upload, give file name]',
                    help='name of tsv with SVid, sample id from vcf, and a CN column')
    parser.add_argument('--collapse-mode',required=True,choices=['gene', 'gene-sv'],metavar='[choose whether collapse to gene level]',
                    help='either collapse rare SV nearby gene, or evaluate as per gene per sv. ')
    parser.add_argument('--remove-control-genes',action=argparse.BooleanOptionalAction,metavar='[--remove-control-genes|--no-remove-control-genes]',
                    help='remove genes with no outliers or not')
    parser.add_argument('--filter-rare',action=argparse.BooleanOptionalAction,metavar='[--filter-rare|--no-filter-rare]',
                    help='remove not rare sv or not')
    parser.add_argument('--flank',type=int,required=True,metavar='[input flank dist]',
                    help='select flank dist to parse')
    parser.add_argument('--zscore-threshold',type=float,required=False,default=3,metavar='[outlier threshold for gene inclusions]',
                    help='select z-score threhsold for selecting genes to include (if we are filtering out genes with no outliers at all. )')
    parser.add_argument('--minimum-support-tissue-count',type=int,required=False,default=5,metavar='[number of outlier signals required to be considered in study]',
                    help='if number of tissues with observed expression for given gene < this number, exclude from study.')

    args = parser.parse_args()

    # read in stuff: 
    # TODO: this part, we want to use the gene_sv_slop.flank.bed file to get info on all SVs involved here. 
    genotypes = pl.scan_csv(args.genotypes,separator='\t',dtypes={'SV': str})
    gene_sv = pl.scan_csv(args.gene_sv,separator='\t',has_header=False,new_columns=['chr','start','end','SV','SVTYPE','Gene'],dtypes={'SV': str})
    gene_sv = gene_sv.with_columns(pl.col('Gene').str.split(by='.').list.first().alias('gene_id'))
    expression_df = pl.scan_csv(args.expressions,separator='\t')
    # extract maf
    if args.maf_mode == 'extract':
        maf_frame = extract_MAF(args.vcf,args.maf_field).rename(columns={args.maf_field:'af'})
    else:
        maf_frame = pd.read_csv(args.maf_file,sep='\t',names=['SV','af'],header=0,dtype={'SV':str})
    maf_frame = pl.from_pandas(maf_frame).lazy()
    # extract length
    if args.length_mode == 'extract':
        length_frame = extract_length(args.vcf,args.length_field).rename(columns={args.length_field:'length'})
    elif args.length_mode == 'upload-SV':
        length_frame = pd.read_csv(args.length_file,sep='\t',names=['SV','length'],header=0,dtype={'SV':str})
    else: 
        length_frame = pd.read_csv(args.length_file,sep='\t',names=['SubjectID','SV','length'],header=0,dtype={'SV':str})

    length_frame = pl.from_pandas(length_frame).lazy()

    # extract copy number
    if args.CN_mode == 'upload':
        CN_frame = pl.scan_csv(args.CN_file,separator='\t',new_columns=['SubjectID','SV','dCN'],dtypes={'SV': str})

    genes = pl.scan_csv(args.genes,separator='\t',has_header=False,new_columns=['chrom','start','end','gene','gene_type','Strand'])
    genes = genes.with_columns(pl.col('gene').str.split(by='.').list.first().alias('gene_id'))

    # process files for getting rare var and expressions
    genotypes = genotypes.filter(pl.col('SVTYPE')!='BND').with_columns(pl.when(pl.col('SVTYPE')=='VNTR').then(pl.lit('DUP')).otherwise(pl.col('SVTYPE')))

    # filter to tissue with at least 5 measurements.
    median_expression_df_protein_coding = expression_df.\
        melt(id_vars=['gene','Ind']).\
        group_by(['gene','Ind']).agg([
            pl.col('value').drop_nans().median().alias('median_outlier'),
            pl.col('value').drop_nans().len().alias('tissue_count')
        ]).rename({args.expression_id_field:'SUBJID','gene':'gene_id'}).\
        with_columns(pl.col('gene_id').str.split(by='.').list.first()).\
        join(genes.select('gene_id'),on='gene_id',how='inner',coalesce=True).collect()
    qualified_expression_protein_coding = median_expression_df_protein_coding.filter(pl.col('tissue_count')>=args.minimum_support_tissue_count).lazy()
    expression_df_protein_coding = qualified_expression_protein_coding.join(expression_df.rename({args.expression_id_field:'SUBJID','gene':'gene_id'}).with_columns([
        pl.col('gene_id').str.split(by='.').list.first()
    ]),how='left',on=['gene_id','SUBJID'],coalesce=True).\
    with_columns(pl.when(pl.col('median_outlier').abs() > args.zscore_threshold).then(1).otherwise(0).alias('Y')).drop(['median_outlier','tissue_count'])
    # next step is to filter to nonzero and high qual alleles.
    # filter to only rare variants or any variants. 
    if args.filter_rare:
        rare_maf_frame = maf_frame.filter((pl.col('af')<0.01)&(pl.col('af')>=0))
    else:
        rare_maf_frame = maf_frame.filter(pl.col('af')>=0)
    true_rare_alleles = genotypes.filter((pl.col("Allele")>0)).join(rare_maf_frame.select("SV"),on='SV',how='inner',coalesce=True)
    # then filter to those that are near by genes,
    noBND_SV_Ind_Gene = true_rare_alleles.join(gene_sv.select(['gene_id','SV']),on='SV',how='inner',coalesce=True).\
        join(expression_df_protein_coding,on=['SUBJID','gene_id'],how='inner',coalesce=True).filter(pl.col('SVTYPE')!='BND')
    # whether to remove control genes or not. 
    if args.remove_control_genes: # continue with groupby tomorrow morning!
        # filter to make sure at least 1 individual is outlier. 
        filtered_total_data = noBND_SV_Ind_Gene.filter(pl.col('Y').sum().over('gene_id')>0)
    else:
        filtered_total_data = noBND_SV_Ind_Gene


    ############################# major differences between normal model and long range #############
    # TODO read in annotation files from all regions in 3 different df lists. 
    # for each annotations, besides the gene, SV column, 
    # we need to rename them with a new suffix: _geneBody, _tssFlank, _tesFlank. 
    # now merge with rest of the annotations:
    genebody_dir = os.path.join(args.annotation_dir,'intermediates')
    tssFlank_dir = os.path.join(args.annotation_dir,'intermediates_tss_flank')
    tesFlank_dir = os.path.join(args.annotation_dir,'intermediates_tes_flank')

    annotation_list_gene_body = glob.glob(os.path.join(genebody_dir,f'*.{args.flank}.tsv'))
    annotation_list_tss_flank = glob.glob(os.path.join(tssFlank_dir,f'*.{args.flank}.tsv'))
    annotation_list_tes_flank = glob.glob(os.path.join(tesFlank_dir,f'*.{args.flank}.tsv'))
    
    # sv annotation by regions:
    svtype_sl = os.path.join(genebody_dir,f'gene_sv_slop.{args.flank}.bed')
    svtype_gb = os.path.join(genebody_dir,f'gene_sv.{args.flank}.bed')
    svtype_sf = os.path.join(tssFlank_dir,f'gene_sv.{args.flank}.bed')
    svtype_ef = os.path.join(tesFlank_dir,f'gene_sv.{args.flank}.bed')
    
    gene_sv_slop_df = pl.read_csv(svtype_sl,
                              separator='\t',
                              has_header=False,
                              new_columns=['chrom','start','end','SV','SVTYPE','Gene']).\
    to_dummies(['SVTYPE'])
    gene_sv_gb_df = pl.read_csv(svtype_gb,
                                separator='\t',
                                has_header=False,
                                new_columns=['chrom','start','end','SV','SVTYPE','Gene']).\
    to_dummies(['SVTYPE'])

    gene_sv_tss_df = pl.read_csv(svtype_sf,
                                separator='\t',
                                has_header=False,
                                new_columns=['chrom','start','end','SV','SVTYPE','Gene']).\
    to_dummies(['SVTYPE'])
    gene_sv_tes_df = pl.read_csv(svtype_ef,
                                separator='\t',
                                has_header=False,
                                new_columns=['chrom','start','end','SV','SVTYPE','Gene']).\
    to_dummies(['SVTYPE'])
    gene_sv_annotation_w_region = gene_sv_slop_df.\
    join(gene_sv_gb_df.select(['Gene','SV',cs.starts_with('SVTYPE_')]),
        on=['Gene','SV'],suffix='_gene_body',how='left',coalesce=True).\
    join(gene_sv_tss_df.select(['Gene','SV',cs.starts_with('SVTYPE_')]),
        on=['Gene','SV'],suffix='_tss_flank',how='left',coalesce=True).\
    join(gene_sv_tes_df.select(['Gene','SV',cs.starts_with('SVTYPE_')]),
        on=['Gene','SV'],suffix='_tes_flank',how='left',coalesce=True).\
    select((~cs.ends_with('_CNV','_DEL','_DUP','_INS','_INV'))).with_columns([
        pl.col('Gene').str.split('.').list.first().alias('GeneName')
    ]).select(pl.exclude(['Gene','chrom','start','end']))
    # combine with other annotations. make sure they have SV gene field. 
    # gene body
    annots_gene_body = []
    for annot_file_gene_body in annotation_list_gene_body:
        annots_gene_body.append(pl.scan_csv(annot_file_gene_body,separator='\t',dtypes={'SV': str}))
    merged_annotations_gene_body = annots_gene_body[0]
    for i in range(1,len(annots_gene_body)):
        merged_annotations_gene_body=merged_annotations_gene_body.join(annots_gene_body[i],how='outer',on=['Gene','SV'],coalesce=True)
    uncollapsed_annotations_gene_body=merged_annotations_gene_body.rename({"Gene":"GeneName"}).\
    with_columns(pl.col('GeneName').str.split('.').list.first().alias('GeneName'))
    # tss flank
    annots_tss_flank = []
    for annot_file_tss_flank in annotation_list_tss_flank:
        annots_tss_flank.append(pl.scan_csv(annot_file_tss_flank,separator='\t',dtypes={'SV': str}))
    merged_annotations_tss_flank = annots_tss_flank[0]
    for i in range(1,len(annots_tss_flank)):
        merged_annotations_tss_flank=merged_annotations_tss_flank.join(annots_tss_flank[i],how='outer',on=['Gene','SV'],coalesce=True)
    uncollapsed_annotations_tss_flank=merged_annotations_tss_flank.rename({"Gene":"GeneName"}).\
    select([
        pl.col(['GeneName','SV']),
        pl.all().exclude(['GeneName','SV']).name.suffix('_tss_flank')
    ]).with_columns(pl.col('GeneName').str.split('.').list.first().alias('GeneName'))

    # tes flank
    annots_tes_flank = []
    for annot_file_tes_flank in annotation_list_tes_flank:
        annots_tes_flank.append(pl.scan_csv(annot_file_tes_flank,separator='\t',dtypes={'SV': str}))
    merged_annotations_tes_flank = annots_tes_flank[0]
    for i in range(1,len(annots_tes_flank)):
        merged_annotations_tes_flank=merged_annotations_tes_flank.join(annots_tes_flank[i],how='outer',on=['Gene','SV'],coalesce=True)
    uncollapsed_annotations_tes_flank=merged_annotations_tes_flank.rename({"Gene":"GeneName"}).\
    select([
        pl.col(['GeneName','SV']),
        pl.all().exclude(['GeneName','SV']).name.suffix('_tes_flank')
    ]).with_columns(pl.col('GeneName').str.split('.').list.first().alias('GeneName'))
    # debug line!!!!!!!!!!!!!
    filtered_total_data.collect().write_csv('test_filter_total_data.tsv',separator='\t',include_header=True)
    # merge all three regions together. 
    uncollapsed_dataset=filtered_total_data.rename({'SUBJID':'SubjectID','gene_id':'GeneName'}).\
        join(uncollapsed_annotations_gene_body,on=['SV','GeneName'],how='left',coalesce=True).\
        join(uncollapsed_annotations_tss_flank,on=['SV','GeneName'],how='left',coalesce=True).\
        join(uncollapsed_annotations_tes_flank,on=['SV','GeneName'],how='left',coalesce=True).\
        join(gene_sv_annotation_w_region.lazy(),on=['SV','GeneName'],how='left',coalesce=True)
    uncollapsed_dataset_with_MAF = uncollapsed_dataset.\
    join(maf_frame.select(['SV','af']),on='SV',how='left',coalesce=True)
    uncollapsed_dataset_with_MAF.collect().write_csv('test_initial_annotations.tsv',separator='\t',include_header=True) # good till here
    # merge in length.
    if args.length_mode == 'upload-SV' or args.length_mode == 'extract':
        length_frame=length_frame.with_columns(pl.col('length').log())
        uncollapsed_dataset_with_MAF_length = uncollapsed_dataset_with_MAF.\
        join(length_frame.select(['SV','length']),on='SV',how='left',coalesce=True)
    else:
        length_frame=length_frame.with_columns(pl.col('length').log())
        uncollapsed_dataset_with_MAF_length = uncollapsed_dataset_with_MAF.\
        join(length_frame.select(['SubjectID','SV','length']),on=['SubjectID','SV'],how='left',coalesce=True)
    # merge in CN. 
    if args.CN_mode == 'upload':
        # replace CN
        uncollapsed_dataset_with_MAF_length_CN = uncollapsed_dataset_with_MAF_length.drop('CN').\
        join(CN_frame.select(['SV','SubjectID','dCN']),on=['SV','SubjectID'],how='left',coalesce=True)
        uncollapsed_dataset_with_MAF_length_CN = uncollapsed_dataset_with_MAF_length_CN.with_columns(pl.col('dCN').fill_null(0))
    else:
        # CN still need special imputations. 
        uncollapsed_dataset_with_MAF_length_CN = uncollapsed_dataset_with_MAF_length
        uncollapsed_dataset_with_MAF_length_CN = uncollapsed_dataset_with_MAF_length_CN.with_columns(
            pl.when((pl.col('CN')==-1)&(pl.col('SVTYPE')=='DUP'))
            .then(pl.col('Allele'))
            .when((pl.col('CN')==-1)&(pl.col('SVTYPE')=='DEL'))
            .then(-pl.col('Allele'))
            .when((pl.col('CN')==-1)&(pl.col('SVTYPE')=='INS'))
            .then(pl.col('Allele'))
            .otherwise(0).alias('dCN')
            )
    uncollapsed_dataset_with_MAF_length_CN = uncollapsed_dataset_with_MAF_length_CN.select(pl.exclude('SVTYPE')).collect()
    print('here!!!')
    # whether to collapse to gene level or not:
    if args.collapse_mode == 'gene-sv':
        uncollapsed_dataset_with_MAF_length_CN=uncollapsed_dataset_with_MAF_length_CN.with_columns(pl.concat_str([pl.col('GeneName'),pl.col('SV')],separator=':').alias('GeneName'))
    # get pvalues:
    uncollapsed_dataset_with_MAF_length_CN=uncollapsed_dataset_with_MAF_length_CN.\
    with_columns(cs.starts_with(args.expression_field).map_alias(lambda c: 'TE_pvalues_' + c.lstrip(args.expression_field)))\
    .select(~cs.starts_with(args.expression_field))
   
    # collapsing annotations. 
    methods=pd.read_csv('collapse_annotation_instructions.tsv',sep='\t')
    collapse_records=methods.to_records()
    collapse_dict={record[1]:record[3] for record in collapse_records}
    del collapse_dict['GeneName']
    del collapse_dict['SubjectID']
    collapse_modified = {}
    for  key,val in collapse_dict.items():
        if key in uncollapsed_dataset_with_MAF_length_CN.columns and key !='dCN':
            collapse_modified[key]=val
        elif key == 'TE_pvalues':
            collapse_modified[key]='first'
        elif key == 'dCN':
            collapse_modified[key]=(lambda x:x.loc[x.abs().idxmax()])


    # iterate through collapse dict as comprehension to do aggregation
    collapsed_complete_dataset = uncollapsed_dataset_with_MAF_length_CN.drop(['SV','Y']).\
        group_by(['SubjectID','GeneName']).agg(
            [
                pl.col(col).max() if method =='max' 
                else pl.col(col).min() if method =='min'
                else pl.when(pl.col(col).is_null().any()).then(pl.col(col).n_unique()-1)
                .otherwise(pl.col(col).n_unique()) if method == 'nunique'
                else cs.starts_with(col).first() if col =='TE_pvalues'
                else pl.col(col).first() if method == 'first'
                else pl.when(pl.col(col).max().abs()>pl.col(col).min().abs())
                .then(pl.col(col).max()).otherwise(pl.col(col).min()) if col == 'dCN'
                else pl.col(col).max()
                for col,method in collapse_modified.items()
            ]
        )
    # leave imputation to after train_test_split
    collapsed_complete_dataset = collapsed_complete_dataset.with_columns(
        cs.starts_with('TE_pvalues').map(lambda x: stats.norm.sf(abs(x))*2*np.sign(x))
    )
    collapsed_complete_dataset.write_csv(args.outfile,separator='\t')

