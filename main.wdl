version 1.0

import "modules/CADD.wdl" as CADD
import "modules/combine_sv_to_gene_roadmaps.wdl" as combine_sv_to_gene_roadmaps
import "modules/extract_rare_variants.wdl" as extract_rare_variants
import "modules/gene_annotations_processing.wdl" as gene_annotations_processing
import "modules/linsight.wdl" as linsight
import "modules/merge_enhancers.wdl" as merge_enhancers
import "modules/phastcon20.wdl" as phastcon20
import "modules/process_roadmaps.wdl" as process_roadmaps
import "modules/sv_to_exon.wdl" as sv_to_exon
import "modules/sv_to_geneABC.wdl" as sv_to_geneABC
import "modules/sv_to_gene_dist.wdl" as sv_to_gene_dist
import "modules/sv_to_gene_enhancers.wdl" as sv_to_gene_enhancers
import "modules/sv_to_gene_flank_processing.wdl" as sv_to_gene_flank_processing
import "modules/sv_to_gene_processing.wdl" as sv_to_gene_processing
import "modules/sv_to_gene_remap.wdl" as sv_to_gene_remap
import "modules/sv_to_gene_tad.wdl" as sv_to_gene_tad
import "modules/combine_annotations_ABC.wdl" as combine_annotations_ABC
import "modules/collect_files.wd" as collect_files
import "modules/vep.wdl" as vep

workflow Watershed_SV {
    input {
        #extract_rare_variants inputs
        File input_vcf
        File metadata
        String pipeline
        Boolean filter_ethnicity
        Boolean filter_rare
    }

    call extract_rare_variants {
        input:
            input_vcf=input_vcf,
            metadata=metadata,
            pipeline=pipeline,
            filter_ethnicity,
            filter_rare
    }

    output {

    }
}
