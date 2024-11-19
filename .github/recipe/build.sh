#!/usr/bin/env bash

set -o xtrace -o nounset -o pipefail -o errexit

mkdir -p ${PREFIX}/libexec/${PKG_NAME}

install_script() {
    script_name=$1
    cp ${SRC_DIR}/scripts/executable_scripts/${script_name} ${PREFIX}/libexec/${PKG_NAME}

tee ${PREFIX}/bin/${script_name} << EOF
    #!/usr/bin/env bash

    exec ${PREFIX}/libexec/${PKG_NAME}/${script_name} "\$@"
EOF
}

export -f install_python_script

python_script_names=(
    combine_all_annotations_ABC_polars.py
    combine_all_annotations_polar.py
    combine_roadmaps.py
    eval_watershed_prep.py
    extract_SV_exon_info.py
    extract_gene_exec.py
    extract_rare_variants.py
    extract_sv_vep_annotations.py
    merge_enhancers.py
    prep_vep_input.py
    sv_to_gene_bw_scores.py
    sv_to_gene_cpg.py
    sv_to_gene_dist.py
    train_test_predict_split_annotation.py
)

shell_script_names=(
    clean_TADS.sh
    extract_sv_vep_annotations.sh
    generate_annotations.sh
    generate_annotations_ABC.sh
    split_bed_by_stateno.sh
    sv_to_gene_roadmap.sh
)

echo ${python_script_names[@]} | tr ' ' '\n' | xargs -I % bash -c 'install_script %'
echo ${shell_script_names[@]} | tr ' ' '\n' | xargs -I % bash -c 'install_script %'

cp ${SRC_DIR}/scripts/executable_scripts/sv_utils.py ${PREFIX}/libexec/${PKG_NAME}
