version 1.0

import "modules/sv_to_gene_bw_scores.wdl" as CADD
import "modules/sv_to_gene_bw_scores.wdl" as GC
import "modules/combine_sv_to_gene_roadmaps.wdl" as combine_sv_to_gene_roadmaps
import "modules/extract_rare_variants.wdl" as extract_rare_variants
import "modules/gene_annotations_processing.wdl" as gene_annotations_processing
import "modules/sv_to_gene_bw_scores.wdl" as linsight
import "modules/merge_enhancers.wdl" as merge_enhancers
import "modules/sv_to_gene_bw_scores.wdl" as phastcon20
import "modules/process_roadmaps.wdl" as process_roadmaps
import "modules/sv_to_exon.wdl" as sv_to_exon
import "modules/sv_to_geneABC.wdl" as sv_to_geneABC
import "modules/sv_to_gene_dist.wdl" as sv_to_gene_dist
import "modules/sv_to_gene_dist.wdl" as sv_to_gene_cpg
import "modules/sv_to_gene_enhancers.wdl" as sv_to_gene_enhancers
import "modules/sv_to_gene_flank_processing.wdl" as sv_to_gene_flank_processing
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

        # gene_annotations_processing
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
        String name_GC
        String stat_method_GC
        Int upper_limit_GC
        Int lower_limit_GC

        # sv_to_gene_bw_scores_linsight
        File bw_linsight
        String name_linsight
        String stat_method_linsight
        Int upper_limit_linsight
        Int lower_limit_linsight

        # sv_to_gene_bw_scores_CADD
        File bw_CADD
        String name_CADD
        String stat_method_CADD
        Int upper_limit_CADD
        Int lower_limit_CADD

        # sv_to_gene_bw_scores_phastcon
        File bw_phastcon
        String name_phastcon
        String stat_method_phastcon
        Int upper_limit_phastcon
        Int lower_limit_phastcon

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

    output {

    }
}
