from pyhpo import Ontology, HPOSet
from pyhpo.stats import EnrichmentModel
import shutil
import os
import ast
import mygene
import pandas as pd
import numpy as np 
def create_gene_names(query):
    genes = []
    for i in query:
        if 'ensembl' in i:
            if type(i['ensembl']) == dict:
                genes.append((i['query'],i['ensembl']['gene']))
            else:
                for j in i['ensembl']:
                    genes.append((i['query'],j['gene']))
    return pd.DataFrame(genes).rename(columns={0:'GeneSymbol',1:'GeneName'})
def hpo_to_genes (alpha,hpo_terms_df):
    _ = Ontology()
    patient_genes = []
    for i,row in hpo_terms_df.iterrows():
        if row['HPO_count_per_UDNid'] != 0:
            # get hpo list
            hpos = ast.literal_eval(row['HPO_ID'])
            # get hpo set
            ci = HPOSet.from_queries(hpos)
            # initiate model
            gene_model = EnrichmentModel('gene')
            # get genes
            genes = gene_model.enrichment(method='hypergeom', hposet=ci)
            genes_for_patient = pd.DataFrame(genes).rename(columns={'enrichment':'enrichment_pval'})
            genes_for_patient['GeneSymbol'] = [i.name for i in genes_for_patient['item']]
            genes_for_patient['Bonferroni_significant_enriched'] = np.where(genes_for_patient.enrichment_pval < alpha/genes_for_patient.shape[0],True,False)
            genes_for_patient['significant_enriched'] = np.where(genes_for_patient.enrichment_pval < alpha,True,False)
            genes_for_patient['UDN_ID'] = row['UDN_ID']
            patient_genes.append(genes_for_patient)
            print(f'{row["UDN_ID"]}')            
    patient_gene_df = pd.concat(patient_genes)
    return patient_gene_df

if __name__ == "__main__":
    hpo_terms = pd.read_csv('/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/Stanford_udn_hpo_ids_sorted_df.csv')
    hpo2gene_df = hpo_to_genes(0.05,hpo_terms)
    mg = mygene.MyGeneInfo()
    genesymbol_list = hpo2gene_df['GeneSymbol'] 
    query_results=mg.querymany(genesymbol_list,scopes='ensembl,symbol', species='human',fields='ensembl.gene,symbol')
    Symbol2Ensembl = create_gene_names(query_results)
    final_hpo2ensembl = hpo2gene_df.merge(Symbol2Ensembl,how='left',on='GeneSymbol')
    final_hpo2ensembl.to_csv('/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/Stanford_udn_hpo_ids_to_geneset.csv',sep='\t',header=True,index=False)


