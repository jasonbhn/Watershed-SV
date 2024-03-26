import pysam
import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='tools to clean up, impute watershed input')
    parser.add_argument('--vcf',type=str,required=True,metavar='[input SV VCF file]',
                    help='select [rare_SV].vcf to parse')
    parser.add_argument('--outfile',type=str,required=True,metavar='[where to output annotation]',
                    help='where to output single indiv vcf. ')
    parser.add_argument('--subjectID',type=str,required=True,metavar='[Sample ID]',
                    help='whose variants are these? ')
    args = parser.parse_args()

    vcf = args.vcf
    outfile = args.outfile
    subjectID = args.subjectID
    unique_var_counter = 0
    with pysam.VariantFile(vcf,'r') as handle:
        with pysam.VariantFile(outfile,'w',header=handle.header) as out:
            out.header.add_meta('INFO', items=[('ID','SubjectID'), 
                                                    ('Number',1), 
                                                    ('Type','String'),
                                                    ('Description','Whose Variant?')])
            handle.header.add_meta('INFO', items=[('ID','SubjectID'), 
                                                    ('Number',1), 
                                                    ('Type','String'),
                                                    ('Description','Whose Variant?')])
            for var in handle.fetch():
                var.info['SubjectID'] = subjectID
                if not var.id:
                    # give unique variant id to unIDed variants. 
                    var.id = f'{var.info["SVTYPE"]}:{str(unique_var_counter)}'
                    unique_var_counter+=1
                if var.info['SVTYPE'] in ['BND','TRA']:
                    # we can't model these yet. 
                    continue
                elif var.info['SVTYPE']=='INV':
                    # due to lack of benchmark, we want to be more confident!
                    if len(var.info['CALLERS'])>=3:
                        out.write(var)
                else:
                    # due to having benchmark, these will be ok. 
                    if len(var.info['CALLERS'])>=2:
                        out.write(var)