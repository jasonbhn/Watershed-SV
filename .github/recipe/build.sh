#!/usr/bin/env bash

set -o xtrace -o nounset -o pipefail -o errexit

PY_VER=$(python3 -V | cut -d ' ' -f 2 | cut -d '.' -f 1,2)
SP_DIR=${PREFIX}/lib/python${PY_VER}/site-packages
mkdir -p ${SP_DIR}
mkdir -p ${PREFIX}/bin

install_script() {
    script_name=$1
    wrapper_name=${script_name//.py/}
    wrapper_name=${wrapper_name//.sh/}
    cp ${SRC_DIR}/scripts/executable_scripts/${script_name} ${PREFIX}/bin/${wrapper_name}
}

export -f install_script

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
    clean_TADs.sh
    generate_annotations.sh
    generate_annotations_ABC.sh
    run_extract_sv_vep_annotations.sh
    split_bed_by_stateno.sh
    sv_to_gene_roadmap.sh
)

echo ${python_script_names[@]} | tr ' ' '\n' | xargs -I % bash -c 'install_script %'
echo ${shell_script_names[@]} | tr ' ' '\n' | xargs -I % bash -c 'install_script %'

cp ${SRC_DIR}/scripts/executable_scripts/sv_utils.py ${SP_DIR}
