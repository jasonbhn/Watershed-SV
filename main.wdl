version 1.0

import "modules/sv_to_gene_bw_scores.wdl" as CADD
import "modules/sv_to_gene_bw_scores.wdl" as GC
import "modules/combine_sv_to_gene_roadmaps.wdl" as combine_sv_to_gene_roadmaps
import "modules/extract_rare_variants.wdl" as extract_rare_variants
import "modules/extract_gene_exec.wdl" as extract_gene_exec
import "modules/sv_to_gene_bw_scores.wdl" as linsight
import "modules/merge_enhancers.wdl" as merge_enhancers
import "modules/sv_to_gene_bw_scores.wdl" as phastcon20
import "modules/process_roadmaps.wdl" as process_roadmaps
import "modules/sv_to_exon.wdl" as sv_to_exon
import "modules/sv_to_geneABC.wdl" as sv_to_geneABC
import "modules/sv_to_gene_dist.wdl" as sv_to_gene_dist
import "modules/sv_to_geneCPG.wdl" as sv_to_gene_cpg
import "modules/sv_to_gene_enhancers.wdl" as sv_to_gene_enhancers
import "modules/sv_to_gene_tss_flank_processing.wdl" as sv_to_gene_flank_tss_processing
import "modules/sv_to_gene_tes_flank_processing.wdl" as sv_to_gene_flank_tes_processing
import "modules/sv_to_gene_slop_processing.wdl" as sv_to_gene_slop_processing
import "modules/sv_to_gene_processing.wdl" as sv_to_gene_processing
import "modules/sv_to_gene_remap.wdl" as sv_to_gene_remap
import "modules/sv_to_gene_tad.wdl" as sv_to_gene_tad
import "modules/combine_annotations_ABC.wdl" as combine_annotations_ABC
import "modules/collect_files.wd" as collect_files
import "modules/vep.wdl" as vep
import "modules/sv_to_geneCPG.wdl" as sv_to_gene_cpg

workflow Watershed_SV {
    input {
        # extract_rare_varians
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
        String TADs_dir

        # merge_enhancers
        File enhancers
        File primary_cells

        # process_roadmaps
        String roadmap_dir

        # vep
        String vep_cache_dir
    }

    call extract_rare_variants {
        input:
            input_vcf=input_vcf,
            metadata=metadata,
            pipeline=pipeline,
            filters=filters,
            filter_ethnicity,
            filter_rare
    }

    call extract_gene_exec {
        input:
            gencode_genes=gencode_genes,
            genome_bound_file=genome_bound_file
    }

    call sv_to_gene_processing {
        input: 
            flank=flank,
            pipeline_bed=extract_rare_variants.pipeline_input,
            genes_bed=extract_genes_exec.genes
    }

    call sv_to_gene_flank_tss_processing {
        input:
            flank=flank,
            genome_bound_file=genome_bound_file,
            genes_bed=extract_genes_exec.genes,
            pipeline_bed=extract_rare_variants.pipeline_input
    }

    call sv_to_gene_flank_tes_processing {
        input:
            flank=flank,
            genome_bound_file=genome_bound_file,
            genes_bed=extract_genes_exec.genes,
            pipeline_bed=extract_rare_variants.pipeline_input
    }

    call sv_to_gene_slop_processing {
        input:
            flank=flank,
            genome_bound_file=genome_bound_file,
            genes_bed=genes_bed,
            pipeline_bed=extract_rare_variants.pipeline_input
    }

    call sv_to_exon {
        input:
            flank=flank,
            exon_bed=extract_gene_exec.exons,
            gene_sv_bed=sv_to_gene_processing.gene_sv_bed
    }

    scatter (sv_file in [sv_to_gene_processing.gene_sv_bed, sv_to_gene_tes_flank_processing.gene_sv_bed, sv_to_gene_tss_flank_processing.gene_sv_bed]) {
        call sv_to_gene_remap {
            input:
                flank=flank,
                remap_crm=remap_crm,
                gene_sv_bed=sv_file
        }

        call sv_to_geneABC {
            input:
                gene_sv_bed=sv_file,
                ABC_enhancers=ABC_enhancers
        }

        call sv_to_gene_cpg {
            input:
                cpg_file=cpg_file,
                flank=flank,
                gene_sv_bed=sv_file
        }
        
        call GC {
            input:
                bw=bw_GC,
                gene_sv_bed=sv_file,
                name="mean_GC_content",
                stat_method="mean",
                upper_limit=100,
                lower_limit=0
        }

        call CADD {
            input:
                bw=bw_CADD,
                gene_sv_bed=sv_file,
                name="top10_CADD",
                stat_method="top10_mean",
                upper_limit=99,
                lower_limit=0
        }

        call linsight {
            input:
                bw=bw_linsight,
                gene_sv_bed=sv_file,
                name="top10_LINSIGHT",
                stat_method="top10_mean",
                upper_limit=1,
                lower_limit=0
        }

        call phastcon20 {
            input:
                bw=bw_phastcon,
                gene_sv_bed=sv_file,
                name="top10_phastCON",
                stat_method="top10_mean",
                upper_limit=1,
                lower_limit=0
        }

        call sv_to_gene_tad {
            input:
                TADs_dir=TADs_dir,
                genome_bound_file=genome_bound_file,
                gene_sv_bed=sv_file
        }

        call merge_enhancers {
            input:
                flank=flank,
                enhancers=enhancers,
                primary_cells=primary_cells,
                gene_sv_bed=sv_file
        }

    }

    call sv_to_gene_dist {
        input:
            flank=flank,
            gene_bed=extract_genes_exec.genes,
            gene_sv_slop_bed=sv_to_gene_slop_processing.sv_gene_slop_bed,
            gene_tss=extract_gene_exec.gene_tss,
            gene_tes=extract_gene_exec.tes
    }

    output {

    }
}
