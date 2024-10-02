#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser(description='Process genomic intervals into single base start and end locations.')
parser.add_argument('intervalBed', metavar='bed-interval', type=str,
                    help='file path to SV genomic intervals')
parser.add_argument('outputBed', metavar='bed-single-base', type=str,
                    help='file path to output decoupled genomic intervals')

args = parser.parse_args()
with open(args.intervalBed,'r') as handle, open(args.outputBed,'w') as out:
    for line in handle:
        chrom,start,end,ID,mapQ = line.rstrip().split('\t')
        start = int(start)
        end = int(end)
        if end - start == 1:
           flag = 'single'
           out.write('\t'.join([chrom,str(start),str(end),ID,mapQ,flag])+'\n')
        else:
           flag_start,flag_end = ('start','end')
           out.write('\t'.join([chrom,str(start),str(start+1),ID,mapQ,flag_start])+'\n')
           out.write('\t'.join([chrom,str(end-1),str(end),ID,mapQ,flag_end])+'\n')


