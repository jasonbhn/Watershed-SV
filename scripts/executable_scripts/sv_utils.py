import pandas as pd
import numpy as np
import pysam
from pysam import VariantFile
import glob
from scipy import stats
# a helper function to apply for filter.
def additional_vcf_info(vcf_path):
    additional_SV_info = {}
    vcf_read = VariantFile(vcf_path,'r')
    for rec in vcf_read.fetch():
        svlen = -1
        af = -1
        if rec.info['SVTYPE']== 'LINE1':
            if rec.info['SVLEN'][0]!= -1:
                svlen = abs(int(rec.info['SVLEN'][0]))
            else:
                svlen = abs(int(rec.info['MEINFO'][2]))
            additional_SV_info[rec.id]= svlen
        elif rec.info['SVTYPE']== 'MEI':
            svlen = abs(int(rec.info['SVLEN'][0]))
            additional_SV_info[rec.id]= svlen
        elif rec.info['SVTYPE']== 'SVA':
            svlen = abs(int(rec.info['SVLEN'][0]))
            additional_SV_info[rec.id]= svlen
        elif rec.info['SVTYPE']== 'ALU':
            svlen = abs(int(rec.info['SVLEN'][0]))
            additional_SV_info[rec.id]= svlen
        elif rec.info['SVTYPE']=='BND':
            svlen=-1
            additional_SV_info[rec.id]= svlen
        elif rec.info['SVTYPE']=='CNV':
            svlen=abs(rec.stop-rec.pos)
            additional_SV_info[f"{rec.id}.DUP_CNV"]= svlen
            additional_SV_info[f"{rec.id}.DEL_CNV"]= svlen 
        elif rec.info['SVTYPE']=='DEL':
            svlen= abs(int(rec.info['SVLEN'][0]))
            additional_SV_info[rec.id]= svlen
        else:
            svlen = abs(rec.info['SVLEN'][0])
            additional_SV_info[rec.id]= svlen
    return additional_SV_info

def extract_race(phenotype_path):
    phenotypes = pd.read_csv(phenotype_path,sep='\t').rename({'SUBJID':'SubjectID'},axis=1).set_index('SubjectID')
    return phenotypes['RACE']
def is_rare(dataframe,MAF_cutoff=0.01):
    AF = dataframe['Allele'].sum()/(dataframe.shape[0]*2)
    if (AF > 0 and AF < MAF_cutoff) or (1-AF >0 and 1-AF <MAF_cutoff):
        return True
    else:
        return False

def maf(dataframe):
    AF = dataframe['Allele'].sum()/(dataframe.shape[0]*2)
    return AF
def svtype_vep(col):
    if col == 'INS':
        return 'INS'
    elif col == 'DEL' or col == 'DEL_CNV':
        return 'DEL'
    elif col == 'DUP':
        return 'DUP'
    elif col == 'DUP_CNV':
        return 'TDUP'
def svtype_annotsv(SVTYPEs):
    mapping= {'DUP_CNV':'CNV', 
              'DEL':'DEL', 
              'CNV':'CNV',
              'DEL_CNV':'CNV', 
              'DUP':'DUP', 
              'ALU':'ALU', 
              'MEI':'DEL', 
              'LINE1':'LINE1',
              'SVA':'SVA',
              'BND':'BND',
              'INV':'INV',
              'INS':'INS',
              'VNTR':'DUP'}
    return [mapping[i] for i in SVTYPEs]
def collapse_types(SVTYPEs):
    mapping= {'DUP_CNV':'DUP_CNV', 
              'DEL':'DEL', 
              'CNV':'CNV',
              'DEL_CNV':'DEL_CNV', 
              'DUP':'DUP',
              'INS':'INS',
              'ALU':'INS', 
              'MEI':'DEL', 
              'LINE1':'INS',
              'SVA':'INS',
              'BND':'BND',
              'INV':'INV',
              'VNTR':'DUP'}
    return [mapping[i] for i in SVTYPEs]
# a helper function to flip the rare variant allele counts...
def is_flipped(dataframe,MAF_cutoff=0.01):
    AF = dataframe['Allele'].sum()/(dataframe.shape[0]*2)
    if 1-AF <MAF_cutoff:
        return True
    else:
        return False
def flipped_type(SVTYPEs):
    flip_dict = {'BND':'BND','DEL_CNV':'DUP_CNV','DUP_CNV':'DEL_CNV','DUP':'DEL','DEL':'INS','INS':'DEL','INV':'INV'}
    return [flip_dict[i] for i in SVTYPEs]
        
    
def no_coord(dataframe,essential_columns=['Gene','SV','Ind','Y']):
    assert pd.Series(essential_columns).isin(dataframe.columns).all(), \
    f"dataframe lack essential columns, it needs to have: {essential_columns}"
    return dataframe[essential_columns]
def vcf_to_numpy_type_lookup(tags,vcf_handle):
    type_list = []
    #vcftype to numpy type:
    vcf_to_npy={'Integer':'int', 
                'Float':'float', 
                'Flag':'object', 
                'Character':'object', 
                'String':'object'}
    # check additional info fields that could be converted. 
    for i in vcf_handle.header.records:
        keys = i.keys()
        if 'ID' in keys and i['ID'] in tags and 'Type' in i.keys():
            type_list.append((i['ID'],vcf_to_npy[i['Type']]))
    tag_list = [i[0] for i in type_list]
    return type_list,tag_list
# helper func
def cn_judge(cn,modecn,gain=True):
    if gain==True:
        if cn > modecn:
            return 2
        else:
            return 0
    else: 
        if cn < modecn:
            return 2
        else:
            return 0
# we assume CNV has already been categorized as rare! now we care about splitting them. 
def extract_genotype_tuples(vcf_path,additional_info_fields=[],
                            filters_pass=['PASS','.'],extract_genotypes=True,
                            filter_ethnicity=False,metadata_path=None):
    # map genotype format to allele
    Alleles_list = []
    SV_bed_list = []
    low_qual_list = []
    gt_dict = {(0,0):0,(0,1):1,(1,1):2,(1,0):1,(0,None):0,(1,None):1,(None,None):-1,(None,):-1,(0,):0,(1,):1}
    # vcf input read with pysam VariantFile class
    vcf_in = VariantFile(vcf_path)
    # get return bedfile numpy matrix dtype
    default_dtypes=[('chrom','object'),('start','int'),
                    ('end','int'),('SV','object'),('SVTYPE','object')]
    # check for dtypes of additional info fields. 
    if len(additional_info_fields)!=0:
        additional_dtypes,sorted_additional_fields=vcf_to_numpy_type_lookup(additional_info_fields,vcf_in)
        dtypes = default_dtypes + additional_dtypes
        
    else:
        dtypes = default_dtypes
    # check if extract genotypes
    if extract_genotypes:
        samples = list(vcf_in.header.samples)
        if filter_ethnicity:
            if metadata_path==None:
                print('please provide metadata containing sample ethnicity')
                return -1
            else:
                race_df = extract_race(metadata_path)
                samples_x = list(set(samples).intersection(set(race_df.index.tolist())))
                samples = list(race_df.loc[samples_x].loc[race_df==3].index)
    for record in vcf_in.fetch():
        record_filters = list(record.filter)
        if len(record_filters) == 0:
            if '.' not in record_filters:
                low_qual_list.append(record.id)
                continue
        if record_filters[0] not in filters_pass:
            low_qual_list.append(record.id)
            continue
        if not record.id:
            continue
        if record.info['SVTYPE']=='TRA':
            continue
        # check for bad chrom id. 
        if 'chr' in record.chrom:
            chrom = record.chrom
        else:
            chrom = 'chr'+record.chrom
        # check if variant is cnv. 
        if record.info['SVTYPE'] != 'CNV':
            if extract_genotypes:
                format_fields = list(record.format.keys())        
                pop_entry = list(map(record.samples.get,samples))
                check_FT = False
                if 'FT' in format_fields:
                    check_FT = True
                check_CN = False
                if 'CN' in format_fields:
                    check_CN = True
                alleles = [gt_dict[i['GT']] if not check_FT else (gt_dict[i['GT']] if i['FT']!='LQ' else -1 )for i in pop_entry]
                CNs = [-1 if not check_CN else int(i['CN']) for i in pop_entry]
                if alleles.count(1)!=0 or alleles.count(2)!=0: 
                    tp=np.array([samples,[record.id]*len(samples),[record.info['SVTYPE']]*len(samples),alleles,CNs]).T
                    Alleles_list.append(tp)
                    if len(additional_info_fields)!=0:
                        SV_bed_list.append(tuple([chrom, record.pos-1, record.stop,
                                        record.id, record.info['SVTYPE']]+\
                                       [record.info[tag] if type(record.info[tag]) is not tuple else record.info[tag][0] for tag in sorted_additional_fields]))
                    else:
                        SV_bed_list.append(tuple([chrom, record.pos-1, record.stop,
                                        record.id, record.info['SVTYPE']]))
            else:
                if len(additional_info_fields)!=0:
                    SV_bed_list.append(tuple([chrom, record.pos-1, record.stop, record.id, record.info['SVTYPE']]+\
                                       [record.info[tag] if type(record.info[tag]) is not tuple else record.info[tag][0] for tag in sorted_additional_fields]))
                else:
                    SV_bed_list.append(tuple([chrom, record.pos-1, record.stop, record.id, record.info['SVTYPE']]))
        else:
            modecn = record.info['MODECN']
            if extract_genotypes:
                format_fields = list(record.format.keys())        
                pop_entry = list(map(record.samples.get,samples))
                check_FT = False
                if 'FT' in format_fields:
                    check_FT = True
                check_CN = False
                if 'CN' in format_fields:
                    check_CN = True
                CNs = [-1 if not check_CN else int(i['CN']) for i in pop_entry]
                # right here we're checking both up and down alleles. But when filtering for rare variants, we should consider both together
                # ie if up + down is not rare, then both should be discarded. 
                up_alleles = [cn_judge(int(i['CN']),modecn) if not check_FT else (cn_judge(int(i['CN']),modecn) if i['FT']!='LQ' else -1) for i in pop_entry]
                down_alleles = [cn_judge(int(i['CN']),modecn,gain=False) if not check_FT else (cn_judge(int(i['CN']),modecn,gain=False) if i['FT']!='LQ' else -1) for i in pop_entry]
                if up_alleles.count(2)!=0:
                    tp=np.array([samples,[record.id+'.DUP_CNV']*len(samples),['DUP_'+record.info['SVTYPE']]*len(samples),up_alleles,CNs]).T
                    Alleles_list.append(tp)
                    if len(additional_info_fields)!=0:
                        SV_bed_list.append(tuple([chrom, record.pos-1, record.stop, record.id+'.DUP_CNV', 'DUP_CNV']+\
                                       [record.info[tag] if type(record.info[tag]) is not tuple else record.info[tag][0] for tag in sorted_additional_fields]))
                    else:
                        SV_bed_list.append(tuple([chrom, record.pos-1, record.stop, record.id+'.DUP_CNV', 'DUP_CNV']))
                if down_alleles.count(2)!=0:
                    tp=np.array([samples,[record.id+'.DEL_CNV']*len(samples),['DEL_'+record.info['SVTYPE']]*len(samples),down_alleles,CNs]).T
                    Alleles_list.append(tp)
                    if len(additional_info_fields)!=0:
                        SV_bed_list.append(tuple([chrom, record.pos-1, record.stop, record.id+'.DEL_CNV', 'DEL_CNV']+\
                                       [record.info[tag] if type(record.info[tag]) is not tuple else record.info[tag][0] for tag in sorted_additional_fields]))
                    else:
                        SV_bed_list.append(tuple([chrom, record.pos-1, record.stop, record.id+'.DEL_CNV', 'DEL_CNV']))
            else:
                if len(additional_info_fields)!=0:
                    SV_bed_list.append(tuple([chrom, record.pos-1, record.stop, record.id, 'CNV']+\
                                   [record.info[tag] if type(record.info[tag]) is not tuple else record.info[tag][0] for tag in sorted_additional_fields]))
                else:
                    SV_bed_list.append(tuple([chrom, record.pos-1, record.stop, record.id, 'CNV']))
    vcf_in.close()
    bed_df=pd.DataFrame(np.array(SV_bed_list,dtype=dtypes))
    return Alleles_list, bed_df
