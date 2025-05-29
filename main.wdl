version 1.0

import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_gene_bw_scores.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/combine_sv_to_gene_roadmaps.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/extract_rare_variants.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/extract_gene_exec.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/merge_enhancers.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/process_roadmaps.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_exon.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_geneABC.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_gene_dist.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_geneCPG.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_gene_enhancers.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_gene_tss_flank_processing.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_gene_tes_flank_processing.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_gene_slop_processing.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_gene_processing.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_gene_remap.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/TAD_annotation_clean_up.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/sv_to_gene_tad.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/combine_annotations_ABC.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/vep_post_processing.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/vep_preprocessing.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/vep_call.wdl"
import "https://raw.githubusercontent.com/jasonbhn/Watershed-SV/refs/heads/WDL/modules/collect_files.wdl"

workflow Watershed_SV {
    input {
        # extract_rare_variants
        File input_vcf
        File metadata
        String pipeline
        Array[String] filters
        Boolean filter_ethnicity
        Boolean filter_rare

        # extract_gene_exec
        File gencode_genes
        File genome_bound_file

        # sv_to_exon
        Int flank

        # sv_to_gene_remap
        File remap_crm

        # sv_to_geneABC
        File ABC_enhancers

        # sv_to_gene_bw_scores_GC
        File bw_GC

        # sv_to_gene_bw_scores_linsight
        File bw_linsight

        # sv_to_gene_bw_scores_CADD
        File bw_CADD

        # sv_to_gene_bw_scores_phastcon
        File bw_phastcon

        # sv_to_gene_cpg
        File cpg_file

        # sv_to_gene_tad
        File TADs_tar

        # merge_enhancers
        File enhancers
        File primary_cells

        # process_roadmaps
        File roadmap_tar

        # vep
        File vep_cache_tar

        # combine_annotations_ABC
        File collapse_instruction
        File gene_list
        File expression_file
        File? maf_file
        File? length_file
        File? CN_file
        Int min_support_tissue
        Float zscore_threshold
        String expression_field
        String expression_id_field
        String maf_mode
        String maf_field
        String length_mode
        String length_field
        String CN_mode
        String collapse_mode
        Boolean remove_control_genes

        # for all modules
        String outdir
        String docker
    }

    call extract_rare_variants.extract_rare_variants {
        input:
            input_vcf=input_vcf,
            metadata=metadata,
            pipeline=pipeline,
            filters=filters,
            filter_ethnicity=filter_ethnicity,
            filter_rare=filter_rare,
            docker=docker,
            memory=32,
            disk_space=40,
            ncpu=4
    }

    call extract_gene_exec.extract_gene_exec {
        input:
            gencode_genes=gencode_genes,
            genome_bound_file=genome_bound_file,
            docker=docker,
            memory=2,
            disk_space=5,
            ncpu=2
    }

    call sv_to_gene_processing.sv_to_gene_processing {
        input: 
            flank=flank,
            pipeline_bed=extract_rare_variants.pipeline_input,
            genes_bed=extract_gene_exec.genes,
            docker=docker,
            memory=8,
            disk_space=10,
            ncpu=4
    }

    call sv_to_gene_tss_flank_processing.sv_to_gene_tss_flank_processing {
        input:
            flank=flank,
            genome_bound_file=genome_bound_file,
            genes_bed=extract_gene_exec.genes,
            pipeline_bed=extract_rare_variants.pipeline_input,
            docker=docker,
            memory=8,
            disk_space=10,
            ncpu=4
    }

    call sv_to_gene_tes_flank_processing.sv_to_gene_tes_flank_processing {
        input:
            flank=flank,
            genome_bound_file=genome_bound_file,
            genes_bed=extract_gene_exec.genes,
            pipeline_bed=extract_rare_variants.pipeline_input,
            docker=docker,
            memory=8,
            disk_space=10,
            ncpu=4
    }

    call sv_to_gene_slop_processing.sv_to_gene_slop_processing {
        input:
            flank=flank,
            genome_bound_file=genome_bound_file,
            genes_bed=extract_gene_exec.genes,
            pipeline_bed=extract_rare_variants.pipeline_input,
            docker=docker,
            memory=8,
            disk_space=10,
            ncpu=4
    }

    call sv_to_exon.sv_to_exon {
        input:
            flank=flank,
            exon_bed=extract_gene_exec.exons,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            docker=docker,
            memory=8,
            disk_space=10,
            ncpu=4
    }

    # remap
    call sv_to_gene_remap.sv_to_gene_remap as sv_to_gene_remap_gene_body {
        input:
            flank=flank,
            remap_crm=remap_crm,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_remap.sv_to_gene_remap as sv_to_gene_remap_tss_flank {
        input:
            flank=flank,
            remap_crm=remap_crm,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_remap.sv_to_gene_remap as sv_to_gene_remap_tes_flank {
        input:
            flank=flank,
            remap_crm=remap_crm,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    # ABC
    call sv_to_geneABC.sv_to_geneABC as sv_to_geneABC_gene_body {
        input:
            flank=flank,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            ABC_enhancers=ABC_enhancers,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_geneABC.sv_to_geneABC as sv_to_geneABC_tss_flank {
        input:
            flank=flank,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            ABC_enhancers=ABC_enhancers,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_geneABC.sv_to_geneABC as sv_to_geneABC_tes_flank {
        input:
            flank=flank,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            ABC_enhancers=ABC_enhancers,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    # cpg
    call sv_to_geneCPG.sv_to_geneCPG as sv_to_geneCPG_gene_body {
        input:
            cpg_file=cpg_file,
            flank=flank,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_geneCPG.sv_to_geneCPG as sv_to_geneCPG_tss_flank {
        input:
            cpg_file=cpg_file,
            flank=flank,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_geneCPG.sv_to_geneCPG as sv_to_geneCPG_tes_flank {
        input:
            cpg_file=cpg_file,
            flank=flank,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    #GC
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as GC_gene_body {
        input:
            flank=flank,
            bw=bw_GC,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            name="mean_GC_content",
            stat_method="mean",
            upper_limit=100,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as GC_tss_flank {
        input:
            flank=flank,
            bw=bw_GC,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            name="mean_GC_content",
            stat_method="mean",
            upper_limit=100,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as GC_tes_flank {
        input:
            flank=flank,
            bw=bw_GC,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            name="mean_GC_content",
            stat_method="mean",
            upper_limit=100,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }

    # CADD
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as CADD_gene_body {
        input:
            flank=flank,
            bw=bw_CADD,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            name="top10_CADD",
            stat_method="top10_mean",
            upper_limit=99,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as CADD_tss_flank {
        input:
            flank=flank,
            bw=bw_CADD,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            name="top10_CADD",
            stat_method="top10_mean",
            upper_limit=99,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as CADD_tes_flank {
        input:
            flank=flank,
            bw=bw_CADD,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            name="top10_CADD",
            stat_method="top10_mean",
            upper_limit=99,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }

    # linsight
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as linsight_gene_body {
        input:
            flank=flank,
            bw=bw_linsight,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            name="top10_LINSIGHT",
            stat_method="top10_mean",
            upper_limit=1,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as linsight_tss_flank {
        input:
            flank=flank,
            bw=bw_linsight,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            name="top10_LINSIGHT",
            stat_method="top10_mean",
            upper_limit=1,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as linsight_tes_flank {
        input:
            flank=flank,
            bw=bw_linsight,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            name="top10_LINSIGHT",
            stat_method="top10_mean",
            upper_limit=1,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    # phastcon20
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as phastcon20_gene_body {
        input:
            flank=flank,
            bw=bw_phastcon,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            name="top10_phastCON",
            stat_method="top10_mean",
            upper_limit=1,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as phastcon20_tss_flank {
        input:
            flank=flank,
            bw=bw_phastcon,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            name="top10_phastCON",
            stat_method="top10_mean",
            upper_limit=1,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_bw_scores.sv_to_gene_bw_scores as phastcon20_tes_flank {
        input:
            flank=flank,
            bw=bw_phastcon,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            name="top10_phastCON",
            stat_method="top10_mean",
            upper_limit=1,
            lower_limit=0,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    # TAD
    call TAD_annotation_clean_up.TAD_annotation_clean_up {
        input:
            TADs_tar=TADs_tar,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_tad.sv_to_gene_tad as sv_to_gene_tad_gene_body {
        input:
            flank=flank,
            cleaned_TAD_beds=TAD_annotation_clean_up.cleaned_TAD_beds,
            genome_bound_file=genome_bound_file,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_tad.sv_to_gene_tad as sv_to_gene_tad_tss_flank {
        input:
            flank=flank,
            cleaned_TAD_beds=TAD_annotation_clean_up.cleaned_TAD_beds,
            genome_bound_file=genome_bound_file,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call sv_to_gene_tad.sv_to_gene_tad as sv_to_gene_tad_tes_flank {
        input:
            flank=flank,
            cleaned_TAD_beds=TAD_annotation_clean_up.cleaned_TAD_beds,
            genome_bound_file=genome_bound_file,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    # enhancers
    call merge_enhancers.merge_enhancers as merge_enhancers_gene_body {
        input:
            flank=flank,
            enhancers=enhancers,
            primary_cells=primary_cells,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call merge_enhancers.merge_enhancers as merge_enhancers_tss_flank {
        input:
            flank=flank,
            enhancers=enhancers,
            primary_cells=primary_cells,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    call merge_enhancers.merge_enhancers as merge_enhancers_tes_flank {
        input:
            flank=flank,
            enhancers=enhancers,
            primary_cells=primary_cells,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=8,
            disk_space=20,
            ncpu=4
    }
    # roadmap
    call process_roadmaps.process_roadmaps as process_roadmaps_gene_body {
        input:
            flank=flank,
            roadmap_tar=roadmap_tar,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=16,
            disk_space=50,
            ncpu=8
    }
    call process_roadmaps.process_roadmaps as process_roadmaps_tss_flank {
        input:
            flank=flank,
            roadmap_tar=roadmap_tar,
            gene_sv_bed=sv_to_gene_tss_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=16,
            disk_space=50,
            ncpu=8
    }
    call process_roadmaps.process_roadmaps as process_roadmaps_tes_flank {
        input:
            flank=flank,
            roadmap_tar=roadmap_tar,
            gene_sv_bed=sv_to_gene_tes_flank_processing.gene_sv_flank_bed,
            docker=docker,
            memory=16,
            disk_space=50,
            ncpu=8
    }


    call sv_to_gene_dist.sv_to_gene_dist {
        input:
            flank=flank,
            gene_bed=extract_gene_exec.genes,
            gene_sv_slop_bed=sv_to_gene_slop_processing.sv_gene_slop_bed,
            gene_tss=extract_gene_exec.gene_tss,
            gene_tes=extract_gene_exec.gene_tes,
            docker=docker,
            memory=8,
            disk_space=5,
            ncpu=4

    }

    call vep_preprocessing.vep_preprocessing {
        input:
            flank=flank,
            gene_sv_slop=sv_to_gene_slop_processing.sv_gene_slop_bed,
            docker=docker,
            memory=16,
            disk_space=32,
            ncpu=8
    }

    call vep_call.vep_call {
        input:
            vep_input=vep_preprocessing.vep_input,
            vep_cache_tar=vep_cache_tar,
            docker="ensemblorg/ensembl-vep:release_109.3",
            memory=16,
            disk_space=300,
            ncpu=8
    }

    call vep_post_processing.vep_post_processing {
        input:
            flank=flank,
            vep_out=vep_call.vep_out,
            docker=docker,
            memory=16,
            disk_space=32,
            ncpu=8
    }
#    call descriptor.collect_files as collect_files_gene_body {
#        input:
#            files=flatten([[sv_to_gene_remap_gene_body.remap_crm_sv_tsv],
#                    [sv_to_geneABC_gene_body.sv_to_gene_flank],
#                    [sv_to_geneCPG_gene_body.sv_to_gene_cpg_dist],
#                    [phastcon20_gene_body.gene_sv_score],
#                    [CADD_gene_body.gene_sv_score],
#                    [GC_gene_body.gene_sv_score],
#                    [linsight_gene_body.gene_sv_score],
#                    [sv_to_gene_tad_gene_body.sv_to_genes_tad],
#                    [merge_enhancers_gene_body.enhancers_by_genes],
#                    process_roadmaps_out_gene_body,
#                    [vep_post_processing.vep],
#                    [sv_to_gene_dist.sv_dist_to_gene],
#                    [sv_to_exon.exon_sv_tsv],
#                    ]),
#            out_tar_basename="gene_body",
#            docker=docker,
#            memory=8,
#            disk_space=60,
#            ncpu=4
#    }
#    call descriptor.collect_files as collect_files_tss_flank {
#        input:
#            files=flatten([[sv_to_gene_remap_tss_flank.remap_crm_sv_tsv],
#                    [sv_to_geneABC_tss_flank.sv_to_gene_flank],
#                    [sv_to_geneCPG_tss_flank.sv_to_gene_cpg_dist],
#                    [phastcon20_tss_flank.gene_sv_score],
#                    [CADD_tss_flank.gene_sv_score],
#                    [GC_tss_flank.gene_sv_score],
#                    [linsight_tss_flank.gene_sv_score],
#                    [sv_to_gene_tad_tss_flank.sv_to_genes_tad],
#                    [merge_enhancers_tss_flank.enhancers_by_genes],
#                    process_roadmaps_out_tss_flank]),
#            out_tar_basename="tss_flank",
#            docker=docker,
#            memory=8,
#            disk_space=60,
#            ncpu=4
#    }
#    call descriptor.collect_files as collect_files_tes_flank {
#        input:
#            files=flatten([[sv_to_gene_remap_tes_flank.remap_crm_sv_tsv],
#                    [sv_to_geneABC_tes_flank.sv_to_gene_flank],
#                    [sv_to_geneCPG_tes_flank.sv_to_gene_cpg_dist],
#                    [phastcon20_tes_flank.gene_sv_score],
#                    [CADD_tes_flank.gene_sv_score],
#                    [GC_tes_flank.gene_sv_score],
#                    [linsight_tes_flank.gene_sv_score],
#                    [sv_to_gene_tad_tes_flank.sv_to_genes_tad],
#                    [merge_enhancers_tes_flank.enhancers_by_genes],
#                    process_roadmaps_out_tes_flank
#                    ]),
#            out_tar_basename="tes_flank",
#            docker=docker,
#            memory=8,
#            disk_space=60,
#            ncpu=4
#    }
    call combine_annotations_ABC.combine_annotations_ABC {
        input:
            collapse_instruction=collapse_instruction,
            flank=flank,
            min_support_tissue=min_support_tissue,
            zscore_threshold=zscore_threshold,
            sv_VCF=input_vcf,
            genotype_VCF=extract_rare_variants.pipeline_input_genotypes,
            gene_list=extract_gene_exec.genes,
            expression_file=expression_file,
            maf_file=maf_file,
            length_file=length_file,
            CN_file=CN_file,
            gene_body_files=flatten([[sv_to_gene_slop_processing.sv_gene_slop_bed],
                    [sv_to_gene_processing.gene_sv_bed],
                    [sv_to_gene_remap_gene_body.remap_crm_sv_tsv],
                    [sv_to_geneABC_gene_body.sv_to_gene_flank],
                    [sv_to_geneCPG_gene_body.sv_to_gene_cpg_dist],
                    [phastcon20_gene_body.gene_sv_score],
                    [CADD_gene_body.gene_sv_score],
                    [GC_gene_body.gene_sv_score],
                    [linsight_gene_body.gene_sv_score],
                    [sv_to_gene_tad_gene_body.sv_to_genes_tad],
                    [merge_enhancers_gene_body.enhancers_by_genes],
                    process_roadmaps_out_gene_body,
                    [vep_post_processing.vep],
                    [sv_to_gene_dist.sv_dist_to_gene],
                    [sv_to_exon.exon_sv_tsv],
                    ]),
            tss_files=flatten([[sv_to_gene_tss_flank_processing.gene_sv_flank_bed],
                    [sv_to_gene_remap_tss_flank.remap_crm_sv_tsv],
                    [sv_to_geneABC_tss_flank.sv_to_gene_flank],
                    [sv_to_geneCPG_tss_flank.sv_to_gene_cpg_dist],
                    [phastcon20_tss_flank.gene_sv_score],
                    [CADD_tss_flank.gene_sv_score],
                    [GC_tss_flank.gene_sv_score],
                    [linsight_tss_flank.gene_sv_score],
                    [sv_to_gene_tad_tss_flank.sv_to_genes_tad],
                    [merge_enhancers_tss_flank.enhancers_by_genes],
                    process_roadmaps_out_tss_flank]),
            tes_files=flatten([[sv_to_gene_tes_flank_processing.gene_sv_flank_bed],
                    [sv_to_gene_remap_tes_flank.remap_crm_sv_tsv],
                    [sv_to_geneABC_tes_flank.sv_to_gene_flank],
                    [sv_to_geneCPG_tes_flank.sv_to_gene_cpg_dist],
                    [phastcon20_tes_flank.gene_sv_score],
                    [CADD_tes_flank.gene_sv_score],
                    [GC_tes_flank.gene_sv_score],
                    [linsight_tes_flank.gene_sv_score],
                    [sv_to_gene_tad_tes_flank.sv_to_genes_tad],
                    [merge_enhancers_tes_flank.enhancers_by_genes],
                    process_roadmaps_out_tes_flank
                    ]),
            expression_field=expression_field,
            expression_id_field=expression_id_field,
            maf_mode=maf_mode,
            maf_field=maf_field,
            length_mode=length_mode,
            length_field=length_field,
            CN_mode=CN_mode,
            collapse_mode=collapse_mode,
            remove_control_genes=remove_control_genes,
            docker=docker,
            memory=64,
            disk_space=60,
            ncpu=32
    }

    output {
        File sv_to_gene_remap_out_gene_body = sv_to_gene_remap_gene_body.remap_crm_sv_tsv
        File sv_to_gene_remap_out_tss_flank = sv_to_gene_remap_tss_flank.remap_crm_sv_tsv
        File sv_to_gene_remap_out_tes_flank = sv_to_gene_remap_tes_flank.remap_crm_sv_tsv
        File sv_to_geneABC_out_gene_body = sv_to_geneABC_gene_body.sv_to_gene_flank
        File sv_to_geneABC_out_tss_flank = sv_to_geneABC_tss_flank.sv_to_gene_flank
        File sv_to_geneABC_out_tes_flank = sv_to_geneABC_tes_flank.sv_to_gene_flank
        File sv_to_geneCPG_out_gene_body = sv_to_geneCPG_gene_body.sv_to_gene_cpg_dist
        File sv_to_geneCPG_out_tss_flank = sv_to_geneCPG_tss_flank.sv_to_gene_cpg_dist
        File sv_to_geneCPG_out_tes_flank = sv_to_geneCPG_tes_flank.sv_to_gene_cpg_dist
        File phastcon20_out_gene_body = phastcon20_gene_body.gene_sv_score
        File phastcon20_out_tss_flank = phastcon20_tss_flank.gene_sv_score
        File phastcon20_out_tes_flank = phastcon20_tes_flank.gene_sv_score
        File CADD_out_gene_body = CADD_gene_body.gene_sv_score
        File CADD_out_tss_flank = CADD_tss_flank.gene_sv_score
        File CADD_out_tes_flank = CADD_tes_flank.gene_sv_score
        File GC_out_gene_body = GC_gene_body.gene_sv_score
        File GC_out_tss_flank = GC_tss_flank.gene_sv_score
        File GC_out_tes_flank = GC_tes_flank.gene_sv_score
        File linsight_out_gene_body = linsight_gene_body.gene_sv_score
        File linsight_out_tss_flank = linsight_tss_flank.gene_sv_score
        File linsight_out_tes_flank = linsight_tes_flank.gene_sv_score
        File sv_to_gene_tad_out_gene_body = sv_to_gene_tad_gene_body.sv_to_genes_tad
        File sv_to_gene_tad_out_tss_flank = sv_to_gene_tad_tss_flank.sv_to_genes_tad
        File sv_to_gene_tad_out_tes_flank = sv_to_gene_tad_tes_flank.sv_to_genes_tad
        File merge_enhancers_out_gene_body = merge_enhancers_gene_body.enhancers_by_genes
        File merge_enhancers_out_tss_flank = merge_enhancers_tss_flank.enhancers_by_genes
        File merge_enhancers_out_tes_flank = merge_enhancers_tes_flank.enhancers_by_genes
        Array[File] process_roadmaps_out_gene_body = process_roadmaps_gene_body.generic_roadmap
        Array[File] process_roadmaps_out_tss_flank = process_roadmaps_tss_flank.generic_roadmap
        Array[File] process_roadmaps_out_tes_flank = process_roadmaps_tes_flank.generic_roadmap
        File combined_annotations = combine_annotations_ABC.combined_dataset
    }
}
