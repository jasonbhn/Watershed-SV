import pandas as pd
import numpy as np
import pyBigWig
import functools
import argparse

def extract_wiggle_scores(var_file,bw_file,score_name,stat_method='max'):
    var = pd.read_csv(var_file,sep='\t',header=None,names=['chrom','start','end','SV','SVTYPE','Gene'])
    score = pyBigWig.open(bw_file)
    chrom_info = score.chroms()
    count = 0
    ids_nan = []
    scores = []
    for i,r in var.iterrows():
        if stat_method!='top10_mean':
           # point max
            if r['chrom'] in chrom_info:
                if r['SVTYPE'] == 'BND' or r['SVTYPE'] == 'INS':
                    stat_score = np.array(score.stats(r['chrom'],r['start']-75,r['end']+75,type=stat_method,exact=True), dtype=np.float64)
                else: 
                    stat_score = np.array(score.stats(r['chrom'],r['start'],r['end'],type=stat_method,exact=True), dtype=np.float64)
                if np.isnan(stat_score):
                    scores.append(-1)
                    ids_nan.append(r['SV'])
                    count+=1
                else:
                    scores.append(stat_score[0])
            else:
                scores.append(np.nan)
        else:
            # average top 10
            if r['chrom'] in chrom_info:
                if r['SVTYPE'] == 'BND' or r['SVTYPE'] == 'INS':
                    interval = np.array(score.values(r['chrom'],r['start']-75,r['end']+75), dtype=np.float64)
                    valid_values = interval[~np.isnan(interval)]
                    if len(valid_values)==0:
                        stat_score = np.nan
                    elif len(valid_values)<=10:
                        stat_score = np.mean(valid_values)
                    else:
                        stat_score = np.mean(valid_values[np.argpartition(valid_values,-10)[-10:]])
                else: 
                    interval = np.array(score.values(r['chrom'],r['start'],r['end']), dtype=np.float64)
                    valid_values = interval[~np.isnan(interval)]
                    if len(valid_values)==0:
                        stat_score = np.nan
                    elif len(valid_values)<=10:
                        stat_score = np.mean(valid_values)
                    else:
                        stat_score = np.mean(valid_values[np.argpartition(valid_values,-10)[-10:]])

                if np.isnan(stat_score):
                    scores.append(-1)
                    ids_nan.append(r['SV'])
                    count+=1
                else:
                    scores.append(stat_score)
            else:
                scores.append(np.nan)
    var[score_name] = scores
    print(f'{score_name}: nan counts = {count}')
    var = var[['SV','Gene',score_name]]
    return var,ids_nan

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gene tp SV distance')
    parser.add_argument('--gene-sv',type=str,required=True,metavar='[gene SV pairing]',
                    help='gene sv pair file')
    parser.add_argument('--in-bigwig',type=str,required=True,metavar='[input for bigwig annotations]',
                    help='path to input bigwig tracks for SV score annotations')
    parser.add_argument('--bigwig-name',type=str,required=True,metavar='[bigwig name]',
                    help='name of input bigwig tracks for SV score annotations')
    parser.add_argument('--stat-method',type=str,required=True,metavar='[method for summarize scores]',
                    help='name of method for summarizing the scores[max,mean,top10_mean]')
    parser.add_argument('--score-upper-limit',type=float,default=None,metavar='[score upper limit]',
                    help='score upper limit, scores > limit = limit (default: %(default)s)') 
    parser.add_argument('--score-lower-limit',type=float,default=None,metavar='[score lower limit]',
                    help='score lower limit, scores < limit = limit (default: %(default)s)') 
    parser.add_argument('--out-gene-sv-score',type=str,required=True,metavar='distance of sv to gene',
                    help='output the sv bigwig file')
    args = parser.parse_args()
    # input argument variables
    gene_sv = args.gene_sv
    in_bigwig = args.in_bigwig
    bigwig_name = args.bigwig_name
    stat_method = args.stat_method
    score_upper_limit = float(args.score_upper_limit)
    score_lower_limit = float(args.score_lower_limit)
    out_gene_sv_score = args.out_gene_sv_score

    # extract scores. 
    score_bigwig,nan_bigwig=extract_wiggle_scores(gene_sv,in_bigwig,bigwig_name,stat_method=stat_method)
    # make sure score values are not over the score limit. 
    score_bigwig.loc[score_bigwig[bigwig_name] > score_upper_limit, bigwig_name] = score_upper_limit
    score_bigwig.loc[score_bigwig[bigwig_name] < score_lower_limit, bigwig_name] = score_lower_limit
    score_bigwig.to_csv(out_gene_sv_score,header=True,index=False,sep='\t')

