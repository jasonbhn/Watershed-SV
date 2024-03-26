import polars as pl
import pandas as pd
import numpy as np
import pyranges as pr
import seaborn as sns
from scipy import stats
from pysam import VariantFile

from matplotlib import pyplot as plt
import statsmodels.formula.api as smf
from itertools import product
import shutil
import os

total_enrichment_matrix = []
for tissue_expression in ['Blood','Fibroblast']:
    for tissue_enhancer in ['Blood','Fibroblast','All']:
        # do enrichment for tissue specific. 
        cohort = "GTEx"
        sv_genotype_old = '/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/new-protein-lincRNA-10k-9.0/intermediates/pipeline_input_genotypes.tsv'
        gene_SV_overlap_cohort_file_old = '/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/new-protein-lincRNA-10k-9.0/intermediates/gene_sv.10000.bed'
        vcf_old = '/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/UDN_MergedSV_vcfs/ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.reIDed.UDN_only.noBNDTRA.vcf.gz'
        expression_file_old = f'/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/UDN_dataset/finalized_UDN_expression/combined_{tissue_expression}_top_60PC_RIN_RUN_SEX_TMM_TPM.nodup.zscore.tsv'
        sv_ABC_overlap_file_old = f'/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/ABC_enhancers/GTEx_SV_{tissue_enhancer}.hg38.region.sorted.bed'
        use_annot_enrichment = False
        annot_col = 'coding_sequence_variant'
        vep_file_old = '/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/new-protein-lincRNA-10k-9.0/intermediates/sv_to_gene_vep.10000.tsv'
        annot_file_old = '/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/new-protein-lincRNA-10k-9.0/intermediates/exon_sv.10000.tsv'

        # stage files to tmp
        sv_genotype = os.path.join('/tmp/',os.path.basename(sv_genotype_old))
        gene_SV_overlap_cohort_file = os.path.join('/tmp/',os.path.basename(gene_SV_overlap_cohort_file_old))
        vcf = os.path.join('/tmp/',os.path.basename(vcf_old))
        expression_file = os.path.join('/tmp/',os.path.basename(expression_file_old))
        sv_ABC_overlap_file = os.path.join('/tmp/',os.path.basename(sv_ABC_overlap_file_old))
        vep_file = os.path.join('/tmp/',os.path.basename(vep_file_old))
        annot_file = os.path.join('/tmp/',os.path.basename(annot_file_old))
        shutil.copyfile(sv_genotype_old,sv_genotype)
        shutil.copyfile(gene_SV_overlap_cohort_file_old,gene_SV_overlap_cohort_file)
        shutil.copyfile(vcf_old,vcf)
        shutil.copyfile(expression_file_old,expression_file)
        shutil.copyfile(sv_ABC_overlap_file_old,sv_ABC_overlap_file)
        shutil.copyfile(vep_file_old,vep_file)
        shutil.copyfile(annot_file_old,annot_file)

        vep_udn = pl.read_csv(vep_file,separator='\t')
        vep_udn = vep_udn.rename({'Gene':'gene'})
        annot_udn = pl.scan_csv(annot_file,separator='\t').fill_null(0)

        UDN_GTEx_blood = pl.scan_csv(expression_file,separator='\t')
        UDN_GTEx_blood = UDN_GTEx_blood.with_columns([pl.when(pl.col('normalized_resid').abs() > 3).then(1).otherwise(0).alias('Y'),
                                                      pl.when(pl.col('Ind').str.contains('GTEX')).then('GTEx').otherwise(cohort).alias('cohort'),
                                                      pl.col('gene').str.split('.').list[0].alias('gene')
                                                     ])
        # UDN global outlier filter
        num_outliers = UDN_GTEx_blood.groupby('Ind').agg(pl.col('Y').sum(),pl.col('cohort').first()).collect()
        num_outliers_cohort_UDN = num_outliers.filter(pl.col('cohort')==cohort)
        num_outliers_cohort_GTEx = num_outliers.filter(pl.col('cohort')=='GTEx')
        IQR_udn = stats.iqr(num_outliers_cohort_UDN['Y'])
        Q3_udn = np.percentile(num_outliers_cohort_UDN['Y'],75)
        Q1_udn = np.percentile(num_outliers_cohort_UDN['Y'],25)
        samples_to_keep = num_outliers_cohort_UDN.filter(pl.col('Y') <= 1.5*IQR_udn+Q3_udn)['Ind']
        UDN_GTEx_blood = UDN_GTEx_blood.filter(pl.col('Ind').is_in(samples_to_keep))
        # read in UDN
        SV_genotypes_cohort=pl.scan_csv(sv_genotype,separator='\t',dtypes={'SV':str})
        gene_SV_overlap_cohort=pl.scan_csv(gene_SV_overlap_cohort_file,separator='\t',dtypes={'SV':str},
                                 new_columns=['chrom','start','end','SV','SVTYPE','gene'])
        gene_SV_overlap_cohort=gene_SV_overlap_cohort.with_columns(pl.col('gene').str.split('.').list.first())
        # read in UDN ABC
        gene_ABC_SV = pl.scan_csv(sv_ABC_overlap_file,
                                        separator='\t',has_header=False).\
                select(['column_6','column_7','column_8','column_4','column_10',
                        'column_11','column_12','column_13',
                        'column_15','column_16']).\
                rename({'column_6':'chrom','column_7':'start','column_8':'end',
                        'column_4':'SV','column_10':'ABC_category',
                        'column_11':'median_ABC.Score','column_12':'max_ABC.Score',
                        'column_13':'sum_ABC.Score','column_15':'percent_celltype','column_16':'gene'})

        cohort_MAF = pl.scan_csv('/oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/Watershed-SV/new-protein-lincRNA-10k-9.0/intermediates/custom_gtex_maf.tsv',
                                 separator='\t',
                                 has_header=True,
                                 dtypes={'SV':str})
        #cohort_MAF = pl.from_pandas(extract_MAF(vcf,fields=['ALNG_AF']).rename(columns={'ALNG_AF':'af'})).lazy()
        # gtex hack lol
        unique_SVs = SV_genotypes_cohort.filter(pl.col('SVTYPE')!='BND').select(pl.col('SV').unique())
        cohort_MAF = unique_SVs.join(cohort_MAF,on='SV',how='left').with_columns(pl.col('af').fill_null(0.5)).lazy()
        rare_svs_cohort=cohort_MAF.with_columns(pl.when(pl.col('af')<0.01).then('rare').otherwise('common').alias('rareness')).select(['SV','af','rareness']).collect()

        rare_genotypes_cohort = SV_genotypes_cohort.filter((pl.col('Allele')>0)).join(rare_svs_cohort.lazy(),on='SV').\
        collect().rename({'SUBJID':'Ind'}).lazy()

        genotype_ind = set(rare_genotypes_cohort.select('Ind').collect()['Ind'])
        expression_ind = set(UDN_GTEx_blood.select('Ind').collect()['Ind'])
        subset_to_use = list(genotype_ind.intersection(expression_ind))

        relevant_outliers=UDN_GTEx_blood.filter(pl.col('Ind').is_in(subset_to_use)).with_columns([pl.when(pl.col('normalized_resid').abs() > i).then(1).otherwise(0).alias(f"Y{i}") for i in range(1,6)])
        relevant_outliers=relevant_outliers.with_columns([pl.when(pl.col(f'Y{i}').sum().over('gene')>0).then(1).otherwise(0).alias(f"Y{i}_background") for i in range(1,6)])

        rare_genotypes_cohort = rare_genotypes_cohort.filter(pl.col('Ind').is_in(subset_to_use))

        coding_annotations = annot_udn.with_columns([
            pl.col('Gene').str.split('.').list[0].alias('gene'),
            pl.when(pl.col('coding_fraction_affected')>0).then(1).otherwise(0).alias('coding_variant')
        ]).\
        select(pl.col(['gene','SV','coding_variant'])).unique()



        rare_sv_enrichment_cohort_ABC = relevant_outliers.\
        join(gene_ABC_SV,on='gene',how='left').\
        join(rare_genotypes_cohort,on=['Ind','SV'],how='left').\
        join(coding_annotations,on=['gene','SV'],how='left').\
        filter(pl.col('cohort')==cohort).fill_null(0).\
        with_columns([
            pl.when(pl.col('rareness')=='rare').then(1).otherwise(0).alias('has_rare_SV'),
            pl.when(pl.col('rareness')=='common').then(1).otherwise(0).alias('has_common_SV'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='DUP')).then(1).otherwise(0).alias('has_rare_DUP'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='DEL')).then(1).otherwise(0).alias('has_rare_DEL'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='INS')).then(1).otherwise(0).alias('has_rare_INS'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='INV')).then(1).otherwise(0).alias('has_rare_INV'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='DUP_CNV')).then(1).otherwise(0).alias('has_rare_DUP_CNV'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='DEL_CNV')).then(1).otherwise(0).alias('has_rare_DEL_CNV'),
                     ]).collect()
        rare_sv_enrichment_cohort_coding = relevant_outliers.\
        join(gene_SV_overlap_cohort,on='gene',how='left').\
        join(rare_genotypes_cohort,on=['Ind','SV'],how='left').\
        join(coding_annotations,on=['gene','SV'],how='left').\
        filter(pl.col('cohort')==cohort).fill_null(0).\
        with_columns([
            pl.when(pl.col('rareness')=='rare').then(1).otherwise(0).alias('has_rare_SV'),
            pl.when(pl.col('rareness')=='common').then(1).otherwise(0).alias('has_common_SV'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='DUP')).then(1).otherwise(0).alias('has_rare_DUP'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='DEL')).then(1).otherwise(0).alias('has_rare_DEL'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='INS')).then(1).otherwise(0).alias('has_rare_INS'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='INV')).then(1).otherwise(0).alias('has_rare_INV'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='DUP_CNV')).then(1).otherwise(0).alias('has_rare_DUP_CNV'),
            pl.when((pl.col('rareness')=='rare')&(pl.col('SVTYPE')=='DEL_CNV')).then(1).otherwise(0).alias('has_rare_DEL_CNV'),
                     ]).collect()
        # we iterate through a few different sum_ABC score threshold.
        varying_ABC_threshold_enrichment_dfs = []
        rare_sv_enrichment_cohort_ABC_variants = rare_sv_enrichment_cohort_ABC.\
        with_columns([
            pl.when((pl.col('has_rare_SV')==True)&(pl.col('coding_variant')==1)).then(1).otherwise(0).alias('has_coding_rare_SV'),
            pl.when((pl.col('has_rare_SV')==True)&
                    (pl.col('coding_variant')==0)&
                    (pl.col('max_ABC.Score')>0)).then(1).otherwise(0).alias('has_ABC_rare_SV'),

        ]).\
        with_columns([
            pl.when((pl.col('has_ABC_rare_SV')==1)&(pl.col('ABC_category')=='genic')).then(1).otherwise(0).alias('has_ABC_genic_rare_SV'),
            pl.when((pl.col('has_ABC_rare_SV')==1)&(pl.col('ABC_category')=='intergenic')).then(1).otherwise(0).alias('has_ABC_intergenic_rare_SV'),
            pl.when((pl.col('has_ABC_rare_SV')==1)&(pl.col('ABC_category')=='promoter')).then(1).otherwise(0).alias('has_ABC_promoter_rare_SV'),
        ]).\
        select([
            "gene","Ind","normalized_resid","cohort","Y1","Y2","Y3","Y4","Y5",
            "Y1_background","Y2_background","Y3_background","Y4_background","Y5_background","Allele",
            "has_rare_SV","has_common_SV","has_rare_DUP","has_rare_DEL","has_rare_INS","has_rare_INV",
            "has_rare_DUP_CNV","has_rare_DEL_CNV",
            "has_ABC_rare_SV","has_coding_rare_SV","has_ABC_genic_rare_SV","has_ABC_intergenic_rare_SV",
            "has_ABC_promoter_rare_SV"
        ])
        rare_sv_enrichment_cohort_other_variants = rare_sv_enrichment_cohort_coding.\
        with_columns([
            pl.when((pl.col('has_rare_SV')==True)&(pl.col('coding_variant')==1)).then(1).otherwise(0).alias('has_coding_rare_SV'),
            pl.lit(0).alias('has_ABC_rare_SV'),
            pl.lit(0).alias('has_ABC_genic_rare_SV'),
            pl.lit(0).alias('has_ABC_intergenic_rare_SV'),
            pl.lit(0).alias('has_ABC_promoter_rare_SV')
        ]).select([
            "gene","Ind","normalized_resid","cohort","Y1","Y2","Y3","Y4","Y5",
            "Y1_background","Y2_background","Y3_background","Y4_background","Y5_background","Allele",
            "has_rare_SV","has_common_SV","has_rare_DUP","has_rare_DEL","has_rare_INS","has_rare_INV",
            "has_rare_DUP_CNV","has_rare_DEL_CNV",
            "has_ABC_rare_SV","has_coding_rare_SV","has_ABC_genic_rare_SV","has_ABC_intergenic_rare_SV",
            "has_ABC_promoter_rare_SV"
        ])

        concat_rare_sv_enrichment_cohort = pl.concat([rare_sv_enrichment_cohort_ABC_variants,rare_sv_enrichment_cohort_other_variants])
        rare_sv_enrichment_cohort_gene = concat_rare_sv_enrichment_cohort.groupby(['gene','Ind']).agg(pl.col("normalized_resid").first(),
                                                                                                      pl.col("cohort").first(),
                                                                                                      pl.col("Y1").first(),
                                                                                                      pl.col("Y2").first(),
                                                                                                      pl.col("Y3").first(),
                                                                                                      pl.col("Y4").first(),
                                                                                                      pl.col("Y5").first(),
                                                                                                      pl.col("Y1_background").first(),
                                                                                                      pl.col("Y2_background").first(),
                                                                                                      pl.col("Y3_background").first(),
                                                                                                      pl.col("Y4_background").first(),
                                                                                                      pl.col("Y5_background").first(),
                                                                                                      pl.col("Allele").max(),
                                                                                                      pl.col("has_rare_SV").max(),
                                                                                                      pl.col("has_common_SV").max(),
                                                                                                      pl.col("has_rare_DUP").max(),
                                                                                                      pl.col("has_rare_DEL").max(),
                                                                                                      pl.col("has_rare_INS").max(),
                                                                                                      pl.col("has_rare_INV").max(),
                                                                                                      pl.col("has_rare_DUP_CNV").max(),
                                                                                                      pl.col("has_rare_DEL_CNV").max(),
                                                                                                      pl.col("has_ABC_rare_SV").max(),
                                                                                                      pl.col("has_coding_rare_SV").max(),
                                                                                                      pl.col("has_ABC_genic_rare_SV").max(),
                                                                                                      pl.col("has_ABC_intergenic_rare_SV").max(),
                                                                                                      pl.col("has_ABC_promoter_rare_SV").max()
                                                                                              )

        enrichment_using_statsmodels = rare_sv_enrichment_cohort_gene.fill_null(0).select(['gene','normalized_resid','Ind',
                                                                                           'Y1','Y2','Y3','Y4','Y5',
                                                                                           'Y1_background','Y2_background',
                                                                                           'Y3_background','Y4_background',
                                                                                           'Y5_background','has_rare_SV','has_common_SV',
                                                                                           'has_rare_DUP','has_rare_DEL','has_rare_INS',
                                                                                           "has_rare_DUP_CNV","has_rare_DEL_CNV",
                                                                                           'has_rare_INV',"has_ABC_rare_SV",
                                                                                           "has_coding_rare_SV","has_ABC_genic_rare_SV",
                                                                                           "has_ABC_intergenic_rare_SV","has_ABC_promoter_rare_SV"
                                                                                           ]).to_pandas()
        enrichment_results = []
        for i in range(1,6):
            # skipping if no outliers. or no SV. 
            skip_condition = enrichment_using_statsmodels.loc[(enrichment_using_statsmodels[f'Y{i}']==1)&
                                                              (enrichment_using_statsmodels['has_ABC_rare_SV']==1)].shape[0]==0
            if skip_condition:
                continue
            model = smf.logit(f'Y{i} ~ has_ABC_rare_SV + has_coding_rare_SV',
                      data = enrichment_using_statsmodels[enrichment_using_statsmodels[f'Y{i}_background']==1]).fit()
            enrichment_results.append(pd.DataFrame({'SVTYPE':['has_ABC_rare_SV',
                                                              'has_coding_rare_SV'],
                                                    'log_OR':[model.params['has_ABC_rare_SV'],
                                                              model.params['has_coding_rare_SV']],
                                                    'SE':[model.bse['has_ABC_rare_SV'],
                                                          model.bse['has_coding_rare_SV']],
                                                    'Z-threshold':i,
                                                    'converged':model.converged,
                                                    'tissue_enhancer':tissue_enhancer,
                                                    'tissue_expression':tissue_expression}))

        if len(enrichment_results)!=0:
            enrichment_df = pd.concat(enrichment_results)
            varying_ABC_threshold_enrichment_dfs.append(enrichment_df)
        print(f'done with tissue: {tissue_enhancer}, {tissue_expression}')
        varying_ABC_threshold_enrichment_df = pd.concat(varying_ABC_threshold_enrichment_dfs)
        total_enrichment_matrix.append(varying_ABC_threshold_enrichment_df)
