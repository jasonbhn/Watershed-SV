from pysam import VariantFile
from sv_utils import *
import pandas as pd
import functools
import os
import numpy as np
import argparse
def impute_dCN (df):
    dup_dCN = np.where((df['CN']==-1)&(df['SVTYPE']=='DUP'),df['Allele'],0)
    del_dCN = np.where((df['CN']==-1)&(df['SVTYPE']=='DEL'),-df['Allele'],0)
    ins_dCN = np.where((df['CN']==-1)&(df['SVTYPE']=='INS'),df['Allele'],0)
    total_dCN = dup_dCN+del_dCN+ins_dCN
    #inversion dcns are automatically 0. 
    return total_dCN

def impute_CN (row):
    if row['CN']==-1:
        if row['DUP']==1:
            return row['Allele']
        elif row['DEL']==1:
            return 2-row['Allele']
        elif row['INS']==1:
            return row['Allele']
        elif row['INV']==1:
            return row['Allele']
    else:
        return row['CN']

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
                    help='directory name of annotations, usually something.../intermediates')
    parser.add_argument('--outfile',type=str,required=True,metavar='[where to output annotation]',
                    help='where to output single set annotation. ')
    parser.add_argument('--expressions',type=str,required=True,metavar='[input gene expression file]',
                    help='select gene expression to parse')
    parser.add_argument('--expression-field',type=str,required=True,metavar='[column name of zscores]',
                    help='select gene expression column to parse')
    parser.add_argument('--expression-id-field',type=str,required=True,metavar='[column name of each sample]',
                    help='select sample column to parse')
    parser.add_argument('--maf-mode',required=True,choices=['upload', 'extract'],metavar='[choose how to get MAF]',
                    help='either upload custom MAF SV->MAF tsv, or extract from info from vcf')
    parser.add_argument('--maf-file',default='None',type=str,metavar='[if --maf-mode upload, give file name]',
                    help='name of tsv with SVid from vcf, and a maf column')
    parser.add_argument('--maf-field',default='Max_AF',type=str,metavar='[if --maf-mode extract, give field name]',
                    help='name of info field from vcf. default follows SVAFotate.')
    parser.add_argument('--length-mode',required=True,choices=['upload', 'extract'],metavar='[choose how to get length]',
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
    parser.add_argument('--flank',type=int,required=True,metavar='[input flank dist]',
                    help='select flank dist to parse')
    args = parser.parse_args()

    # read in stuff
    genotypes = pd.read_csv(args.genotypes,sep='\t',dtype={'SV': str})
    gene_sv = pd.read_csv(args.gene_sv,sep='\t',header=None,names=['chr','start','end','SV','SVTYPE','Gene'],dtype={'SV': str})
    gene_sv['gene_id'] = gene_sv['Gene'].str.split('.').str[0] # immutable with respect to ensembl version
    expression_df = pd.read_csv(args.expressions,sep='\t')

    # extract maf
    if args.maf_mode == 'extract':
        maf_frame = extract_MAF(args.vcf,args.maf_field).rename(columns={args.maf_field:'af'})
    else:
        maf_frame = pd.read_csv(args.maf_file,sep='\t',names=['SV','af'],header=0)
    maf_frame['SV'] = maf_frame.SV.astype('str')
    # extract length
    if args.length_mode == 'extract':
        length_frame = extract_length(args.vcf,args.length_field).rename(columns={args.length_field:'length'})
    else:
        length_frame = pd.read_csv(args.length_file,sep='\t',names=['SV','length'],header=0)
    length_frame['SV'] = length_frame.SV.astype('str')
    # extract copy number
    if args.CN_mode == 'upload':
        CN_frame = pd.read_csv(args.CN_file,sep='\t',names=['SubjectID','SV','dCN'],header=0)
        CN_frame['SV'] = CN_frame.SV.astype('str')

    genes = pd.read_csv(args.genes,sep='\t',header=None,names=['chrom','start','end','gene','gene_type','Strand'])
    genes['gene_id'] = genes['gene'].str.split('.').str[0] # immutable with respect to ensembl version
    # process files for getting rare var and expressions
    expression_df['Y'] = np.where(expression_df[args.expression_field].abs() > 3.0,'outlier','control')
    genotypes = genotypes[genotypes['SVTYPE']!='BND']
    expression_df['gene_id'] = expression_df['gene'].str.split('.').str[0]
    expression_df_protein_coding = expression_df[expression_df['gene_id'].isin(genes.gene_id)]
    # next step is to filter to nonzero and high qual alleles.
    rare_maf_frame = maf_frame[(maf_frame.af <0.01)&(maf_frame.af>=0)]
    true_rare_var = genotypes[genotypes.SV.isin(rare_maf_frame.SV)]
    true_rare_alleles = true_rare_var[(true_rare_var['Allele']>0)].copy()
    expression_df_protein_coding.rename(columns={'gene':'Gene'},inplace=True)

    # then filter to those that are near by genes,
    SV_Ind_Gene_noExpn = true_rare_alleles.merge(gene_sv[['gene_id','SV']],on=['SV'],how='inner')
    # SV->Ind->Gene(further see which one has measured expression)
    SV_Ind_Gene = SV_Ind_Gene_noExpn.merge(expression_df_protein_coding.rename(columns={args.expression_id_field:'SUBJID'}),on=['SUBJID','gene_id'],how='inner')


    # filter to no breakend or not? currently let's just remove bnd. 
    cond_noBND=(SV_Ind_Gene['SVTYPE'] !='BND')
    noBND_SV_Ind_Gene = SV_Ind_Gene[cond_noBND]
    if args.remove_control_genes:
        # filter to make sure at least 1 individual is outlier. 
        gene_keep_condition_by_Y = noBND_SV_Ind_Gene.groupby('Gene').apply(lambda x: True if np.count_nonzero(x['Y'].to_numpy()=='outlier') >0 else False)
        gene_kept = gene_keep_condition_by_Y[gene_keep_condition_by_Y==True].index.unique()
        index_kept = noBND_SV_Ind_Gene.Gene.isin(gene_kept)
        filtered_total_data = noBND_SV_Ind_Gene.loc[index_kept,:]
    else:
        filtered_total_data = noBND_SV_Ind_Gene

    # now merge with rest of the annotations:
    annotation_list = glob.glob(os.path.join(args.annotation_dir,f'*.{args.flank}.tsv'))
    # combine with other annotations. make sure they have SV gene field. 
    annots = []
    for annot_file in annotation_list:
        annots.append(pd.read_csv(annot_file,sep='\t',dtype={'SV': str}))
    merged_annotations=functools.reduce(lambda left,right: pd.merge(left,right,how='outer',on=['Gene','SV']),annots)
    uncollapsed_annotations=merged_annotations.rename(columns={"Gene":"GeneName","Ind":"SubjectID"})
    uncollapsed_annotations['GeneName']=uncollapsed_annotations['GeneName'].str.split('.').str[0]
    #uncollapsed_annotations.rename(columns={'SV':'SVid'},inplace=True)
    #,'SV':'SVid'}).\
    uncollapsed_dataset=filtered_total_data.rename(columns={'SUBJID':'SubjectID','gene_id':'GeneName'}).\
        merge(uncollapsed_annotations,on=['SV','GeneName'],how='left')
    # merge in maf.  #rename(columns={'SVid':'SV'}).\
    uncollapsed_dataset_with_MAF = uncollapsed_dataset.\
    merge(maf_frame[['SV','af']],on='SV',how='left')
    # merge in length.
    length_frame['length'] = np.log(length_frame['length'])
    uncollapsed_dataset_with_MAF_length = uncollapsed_dataset_with_MAF.\
    merge(length_frame[['SV','length']],on='SV',how='left')
    # merge in CN. 
    if args.CN_mode == 'upload':
        # replace CN
        uncollapsed_dataset_with_MAF_length_CN = uncollapsed_dataset_with_MAF_length.drop('CN',axis=1).\
        merge(CN_frame[['SV','SubjectID','dCN']],on=['SV','SubjectID'],how='left')
        uncollapsed_dataset_with_MAF_length_CN['dCN'] = uncollapsed_dataset_with_MAF_length_CN['dCN'].fillna(0)
    else:
        # CN still need special imputations. 
        uncollapsed_dataset_with_MAF_length_CN = uncollapsed_dataset_with_MAF_length
        uncollapsed_dataset_with_MAF_length_CN['dCN'] = impute_dCN(uncollapsed_dataset_with_MAF_length_CN)
    # dummify some stuff.
    uncollapsed_dataset_with_MAF_length_CN_SVTYPE = pd.get_dummies(uncollapsed_dataset_with_MAF_length_CN,columns=['SVTYPE'])
    # do some cleaning for gene_id
    uncollapsed_dataset_with_MAF_length_CN = uncollapsed_dataset_with_MAF_length_CN_SVTYPE.drop('Gene',axis=1)
    # imputations. 
    methods=pd.read_csv('collapse_annotation_instructions.tsv',sep='\t')
    collapse_records=methods.to_records()
    collapse_dict={record[1]:record[3] for record in collapse_records}
    del collapse_dict['GeneName']
    del collapse_dict['SubjectID']
    
    # whether to collapse to gene level or not:
    if args.collapse_mode == 'gene-sv':
        uncollapsed_dataset_with_MAF_length_CN['GeneName']=uncollapsed_dataset_with_MAF_length_CN['GeneName']+':'+uncollapsed_dataset_with_MAF_length_CN['SV']
    # get pvalues:
    uncollapsed_dataset_with_MAF_length_CN=uncollapsed_dataset_with_MAF_length_CN.rename(columns={args.expression_field:'TE_pvalues'})
    collapse_modified = {}
    for  key,val in collapse_dict.items():
        if key in uncollapsed_dataset_with_MAF_length_CN.columns and key !='dCN':
            collapse_modified[key]=val
        elif key == 'dCN':
            collapse_modified[key]=(lambda x:x.loc[x.abs().idxmax()])
    
    collapsed_complete_dataset=uncollapsed_dataset_with_MAF_length_CN.drop(['SV','Y'],axis=1).\
        groupby(['SubjectID','GeneName']).\
            aggregate(collapse_modified) # in this way we correctly impute dCN
    collapsed_complete_dataset.reset_index(inplace=True)
    # leave imputation to after train_test_split
    collapsed_complete_dataset['TE_pvalues'] = stats.norm.sf(abs(collapsed_complete_dataset['TE_pvalues']))*2*np.sign(collapsed_complete_dataset['TE_pvalues'])
    collapsed_complete_dataset.to_csv(args.outfile,sep='\t',header=True,index=False)

