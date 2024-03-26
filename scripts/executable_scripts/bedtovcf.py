from vcftobed import BED_DICT,BED,VCF
import argparse


if __name__ == "__main__":
    # get args
    parser = argparse.ArgumentParser(description='VCF modification based on (post-liftover)BED')
    parser.add_argument('--vcf-in',type=str,metavar='[input VCF file]',
                        help='select .vcf to modify')
    parser.add_argument('--bed',type=str,metavar='[input bed file]',
                        help='select bed file for modification')
    parser.add_argument('--vcf-out',type=str,metavar='[output vcf file]',
                        help='specify vcf out file path')
    args = parser.parse_args()
    # file paths
    vcfin_path=args.vcf_in
    vcfout_path=args.vcf_out
    bed_path=args.bed
    # initialize bed-dict to take note of all vcf entries with successful liftover
    bed_dict = BED_DICT(bed_path)
    out = []
    with open(vcfin_path,'r') as handle:
        for line in handle:
            if line[0]=='#':
                out.append(line.rstrip())
            else:
                vcf = VCF(line)
                if vcf.modify_from_bed(bed_dict):
                    out.append(str(vcf))
    with open(vcfout_path,'w') as handle:
        for line in out:  
            handle.write(line+'\n')


