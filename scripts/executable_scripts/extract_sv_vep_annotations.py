import pandas as pd
import numpy as np
import sys
tmp = sys.argv[1]
vep_out = sys.argv[2]
VEP_out=pd.read_csv(tmp,
                    sep='\t',
                    comment='#',
                    names=['VarID','gene','Feature_Type','Consequence','IMPACT'])
#tmp1=pd.concat([VEP_out, VEP_out.Consequence.str.get_dummies(sep=',')],axis=1).drop('Consequence',axis=1)
#tmp2=pd.concat([tmp1,tmp1.IMPACT.str.get_dummies(sep=',')],axis=1).drop('IMPACT',axis=1)
#tmp3=tmp2.drop(['Feature_Type'],axis=1)

#vep_processed_annotations = tmp3.groupby('VarID').aggregate('max').reset_index()
#vep_processed_annotations['SV'] = vep_processed_annotations.VarID.str.split(':').str[:-1].str.join(':')
#vep_processed_annotations['Gene'] = vep_processed_annotations.VarID.str.split(':').str[-1]
#vep_processed_annotations['SV'] = vep_processed_annotations.SV.astype('str')
tmp1=pd.concat([VEP_out, VEP_out.Consequence.str.get_dummies(sep=',')],axis=1).drop('Consequence',axis=1)
tmp2=pd.concat([tmp1,tmp1.IMPACT.str.get_dummies(sep=',')],axis=1).drop('IMPACT',axis=1)

tmp2['Gene'] = tmp2.VarID.str.split(':').str[-1].str.split('.').str[0]
tmp2_gene = tmp2[(tmp2.gene!='-')&(tmp2.gene==tmp2.Gene)]
tmp2_regulatory = tmp2[tmp2.gene=='-']
vep_processed_gene_annotations = tmp2_gene.groupby('VarID').aggregate('max').reset_index()
vep_processed_regulatory_annotations = tmp2_regulatory.groupby('VarID').aggregate('max').reset_index()
vep_processed_regulatory_annotations['gene'] = vep_processed_regulatory_annotations['Gene']
# coding and intron doesn't make sense.
vep_processed_gene_annotations['intron_variant']= np.where(vep_processed_gene_annotations['coding_sequence_variant']==1,
                                                           0,
                                                           vep_processed_gene_annotations['intron_variant'])
vep_processed_annotations = pd.concat([vep_processed_gene_annotations,vep_processed_regulatory_annotations])
vep_processed_annotations['SV'] = vep_processed_annotations.VarID.str.split(':').str[:-1].str.join(':')
vep_processed_annotations['SV'] = vep_processed_annotations.SV.astype('str')
vep_processed_annotations_final = vep_processed_annotations.groupby(['VarID','Gene']).aggregate('max').reset_index().drop(['VarID','Feature_Type','gene'],axis=1)

vep_processed_annotations_final.to_csv(vep_out,sep='\t',header=True,index=False)
