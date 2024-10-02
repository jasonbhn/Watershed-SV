#!/usr/bin/env python3

import pysam
import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='tools to clean up, impute watershed input')
    parser.add_argument('--vcf',type=str,required=True,metavar='[input SV VCF file]',
                    help='select [rare_SV].vcf to parse')
    parser.add_argument('--outfile',type=str,required=True,metavar='[where to output annotation]',
                    help='where to output single indiv vcf. ')
    
    args = parser.parse_args()

    vcf = args.vcf
    outfile = args.outfile
    with pysam.VariantFile(vcf,'r') as handle:
        with pysam.VariantFile(outfile,'w',header=handle.header) as out:
            out.header.add_meta('INFO', items=[('ID','SEQ'), 
                                                    ('Number',1), 
                                                    ('Type','String'),
                                                    ('Description','Sequence for insertion')])
            handle.header.add_meta('INFO', items=[('ID','SEQ'), 
                                                    ('Number',1), 
                                                    ('Type','String'),
                                                    ('Description','Sequence for insertion')])
            for var in handle.fetch():
                if var.chrom in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']:
                    if var.info['SVTYPE']=='INS':
                        sample = var.info['SubjectID']
                        var.info.clear()
                        if var.alts[0]!='<INS>':
                            var.info['SVTYPE'] = 'INS'
                            var.info['SubjectID'] = sample
                            var.info['SEQ'] = var.alts[0]
                            out.write(var)
                    elif var.info['SVTYPE']=='INV':
                        sample = var.info['SubjectID']
                        var.info.clear()
                        if var.alts[0]!='<INV>':
                            var.alts=('<INV>',)
                        var.info['SVTYPE']='INV'
                        var.info['SubjectID'] = sample
                        out.write(var)

                    elif var.info['SVTYPE']=='DUP':
                        sample = var.info['SubjectID']
                        var.info.clear()
                        if var.alts[0]!='<DUP>':
                            var.alts=('<DUP>',)
                        var.info['SVTYPE']='DUP'
                        var.info['SubjectID'] = sample
                        out.write(var)
                    elif var.info['SVTYPE']=='DEL':
                        sample = var.info['SubjectID']
                        var.info.clear()
                        if var.alts[0]!='<DEL>':
                            var.alts=('<DEL>',)
                        var.info['SVTYPE']='DEL'
                        var.info['SubjectID'] = sample
                        out.write(var)
