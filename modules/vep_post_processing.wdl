
version 1.0

task vep_post_processing {
    input{
        Int flank
        File vep_out
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output{
        File vep = "sv_to_gene_vep.~{flank}.tsv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
    command <<<
        extract_sv_vep_annotations "~{vep_out}" "sv_to_gene_vep.~{flank}.tsv"
    >>>
}
