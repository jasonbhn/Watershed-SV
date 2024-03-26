from pysam import VariantFile

def extract_bed(vcf_path,bed_path):
    vcf_in = VariantFile(vcf_path)   
    with open(bed_path,'w') as handle:
        for record in vcf_in.fetch():
            if list(record.filter)[0]!='PASS' or record.info['SVTYPE']=='BND':
                continue            
            handle.write(f'{record.contig}\t{record.pos-1}\t{record.stop}\t{record.id}\t{record.info["SVTYPE"]}\t{abs(record.info["SVLEN"][0])}\t{record.info["MAF"]}\n') 

if __name__ == "__main__":
    # get args
    parser = argparse.ArgumentParser(description='VCF ')
    parser.add_argument('--vcf',type=str,metavar='[input VCF file]',
                        help='select .vcf to convert')
    parser.add_argument('--bed',type=str,metavar='[output bed file]',
                        help='convert to which bed file')
    args = parser.parse_args()
    # file paths
    vcf_path=args.vcf
    bed_path=args.bed
    extract_bed(vcf_path,bed_path)
