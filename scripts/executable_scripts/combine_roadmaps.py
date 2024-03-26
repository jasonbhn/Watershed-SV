import pandas as pd
import functools
import io
import os
import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='merge enhancer annotations')
    parser.add_argument('--gene-sv-roadmap-dir',type=str,required=True,metavar='[roadmap state-specific gene sv directory]',
                    help='roadmap state-specific gene sv directory')
    parser.add_argument('--out-combined-roadmap',type=str,required=True,metavar='[out combined roadmap annotations]',
                    help='out_combined_roadmap')
    args = parser.parse_args()
    # input argument variables
    gene_sv_roadmap_dir = args.gene_sv_roadmap_dir
    out_combined_roadmap = args.out_combined_roadmap

    string="""
    STATE NO.	MNEMONIC	DESCRIPTION	COLOR NAME	COLOR CODE
    1	TssA	Active TSS	Red	255,0,0
    2	PromU	Promoter Upstream TSS	Orange Red	255,69,0
    3	PromD1	Promoter Downstream TSS 1	Orange Red	255,69,0
    4	PromD2	Promoter Downstream TSS 2	Orange Red	255,69,0
    5	Tx5	Transcribed - 5' preferential	Green	0,128,0
    6	Tx	Strong transcription	Green	0,128,0
    7	Tx3	Transcribed - 3' preferential	Green	0,128,0
    8	TxWk	Weak transcription	Lighter Green	0,150,0
    9	TxReg	Transcribed & regulatory (Prom/Enh)	Electric Lime	194,225,5
    10	TxEnh5	Transcribed 5' preferential and Enh	Electric Lime	194,225,5
    11	TxEnh3	Transcribed 3' preferential and Enh	Electric Lime	194,225,5
    12	TxEnhW	Transcribed and Weak Enhancer	Electric Lime	194,225,5
    13	EnhA1	Active Enhancer 1	Orange	255,195,77
    14	EnhA2	Active Enhancer 2	Orange	255,195,77
    15	EnhAF	Active Enhancer Flank	Orange	255,195,77
    16	EnhW1	Weak Enhancer 1	Yellow	255,255,0
    17	EnhW2	Weak Enhancer 2	Yellow	255,255,0
    18	EnhAc	Primary H3K27ac possible Enhancer	Yellow	255,255,0
    19	DNase	Primary DNase	Lemon	255,255,102
    20	ZNF_Rpts	ZNF genes & repeats	Aquamarine	102,205,170
    21	Het	Heterochromatin	Light Purple	138,145,208
    22	PromP	Poised Promoter	Pink	230,184,183
    23	PromBiv	Bivalent Promoter	Dark Purple	112,48,160
    24	ReprPC	Repressed Polycomb	Gray	128,128,128
    25	Quies	Quiescent/Low	White	255,255,255
    """
    states_25_detail=pd.read_csv(io.StringIO(string),sep='\t')
    allstates=[]
    states_num=[]
    for i in range(1,26):
        roadmap_state_file = os.path.join(gene_sv_roadmap_dir,f'roadmap_multitissue_sv_to_gene.{i}.tsv')
        if os.path.exists(roadmap_state_file):
            # only those that exists. 
            allstates.append(pd.read_csv(roadmap_state_file,sep='\t',dtype={'SV':str,'Gene':str}))
            # now we record the state we extracted. 
            states_num.append(i-1)

    allstate_df=functools.reduce(lambda left,right: pd.merge(left,right,how='outer',on=['SV','Gene']),allstates).fillna(0)

    allstate_df.columns=['SV','Gene']+[states_25_detail.MNEMONIC.tolist()[i] for i in states_num]

    allstate_df.to_csv(out_combined_roadmap,sep='\t',header=True,index=False)
