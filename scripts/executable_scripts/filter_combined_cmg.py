#!/usr/bin/env python3

import argparse
from pysam import VariantFile
def preproc_paragraph(vcf,out_vcf):
    with VariantFile(vcf) as handle:
        handle.header.add_meta('INFO', items=[('ID','MEICLASS'),
                                            ('Number',1),
                                            ('Type','String'),
                                            ('Description','MEI subclass')])
        with VariantFile(out_vcf,'w',header=handle.header) as out:
            for record in handle.fetch():
                out_record = record.copy()
                if type(record.info['SVTYPE'])!=str:
                    return 1
                if record.info['SVTYPE']=="MEI":
                    out_record.info['SVTYPE']="DEL"
                    out_record.alleles=("N","<DEL>")
                    out_record.ref="N"
                    out_record.alts = ('<DEL>',)
                    out_record.info['MEICLASS']=out_record.alts[0]
                elif record.info['SVTYPE']!='INS':
                    #make sure the ref alt is N and <TYPE>
                    out_record.alleles=('N',f'<{out_record.info["SVTYPE"]}>')
                    out_record.ref = 'N'
                    out_record.alts = (f'<{out_record.info["SVTYPE"]}>',)
                out.write(out_record)
                

                    
def filter_cmg_combined(vcf,out_vcf):
    with VariantFile(vcf) as handle:
        with VariantFile(out_vcf,'w',header=handle.header) as out:
            out.header.add_meta('INFO', items=[('ID','CHR2'),
                                                ('Number',1),
                                                ('Type','String'),
                                                ('Description','Complex SV second chromosome')])
            out.header.add_meta('INFO', items=[('ID','CALLERS'),
                                                ('Number','.'),
                                                ('Type','String'),
                                                ('Description','CALLERS consensus')])
            out.header.add_meta('FORMAT', items=[('ID','SP'),
                                                ('Number','.'),
                                                ('Type','String'),
                                                ('Description','CALLERS')])
            consensus=0
            svtools=0
            multiparliament2=0
            for record in handle.fetch():
                if len(record.info['IDLIST'])>1:
                    consensus+=1
                    out.write(record)
                elif 'CALLERS' not in record.info:
                    svtools+=1
                    out.write(record)
                elif 'CALLERS' in record.info:
                    if type(record.info['CALLERS'])==tuple:
                        multiparliament2+=1
                        out.write(record)
                    elif record.info['SVTYPE']=='INS' and type(record.info['CALLERS'])==str:
                        if record.info['CALLERS']=='MANTA':
                            multiparliament2+=1
                            out.write(record)
            print(f'{consensus} in consensus\n{svtools} in svtools only\n{multiparliament2} supported by multiple callers in parliament2')
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF file from combined CMG individual")
    parser.add_argument("--input_vcf", required=True, type=str, help="Input VCF file name")
    parser.add_argument("--output_vcf", required=True, type=str, help="Output VCF file name")
    args = parser.parse_args()

    input_vcf_filename = args.input_vcf
    output_vcf_filename = args.output_vcf
    preproc_paragraph(input_vcf_filename, output_vcf_filename)
