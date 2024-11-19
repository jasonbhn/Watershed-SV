#!/usr/bin/env bash

set -o xtrace -o nounset -o pipefail -o errexit

mkdir -p ${PREFIX}/libexec/${PKG_NAME}

install_script() {
    script_name=$1
    cp ${SRC_DIR}/scripts/executable_scripts/${script_name}.py ${PREFIX}/libexec/${PKG_NAME}
    chmod +x ${PREFIX}/libexec/${PKG_NAME}/${script_name}.py 
    sed -i '1i #!/usr/bin/env python3' ${PREFIX}/libexec/${PKG_NAME}/${script_name}.py

tee ${PREFIX}/bin/${script_name} << EOF
    #!/usr/bin/env bash

    exec ${PREFIX}/libexec/${PKG_NAME}/${script_name}.py "\$@"
EOF
}

export -f install_script

script_names=(
    #bedtovcf  missing vcftobed - is it supposed to be vcftobed_v2?
    collapse_annotation
    combine_all_annotations
    combine_all_annotations_ABC_polars
    combine_all_annotations_polar
    combine_raw_annotations
    combine_roadmaps
    convert_vcf_to_paragraph_input
    enrichment_by_rare_SV
    eval_watershed_prep
    #extract_SV_annotations  missing some arguments
    extract_SV_exon_info
    extract_gene_exec
    extract_rare_variants
    filterByCaller
    filter_combined_cmg
    filter_splice_clusters
    interval_to_pair
    leafcutter_cluster_regtools
    merge_enhancers
    #outlier_calling_exec  has typo on line 12 - parser.add_arl
    regress_out_PC
    regress_out_PC_no_age
    relevant_gene_enrichment_for_outliers
    sv_to_gene_bw_scores
    sv_to_gene_cpg
    sv_to_gene_dist
    train_test_predict_split_annotation
    #vcftobed_v2  missing import argparse
)

echo ${script_names[@]} | tr ' ' '\n' | xargs -I % bash -c 'install_script %'

cp ${SRC_DIR}/scripts/executable_scripts/sv_utils.py ${PREFIX}/libexec/${PKG_NAME}
