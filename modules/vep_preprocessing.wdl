
version 1.0

task vep_preprocessing {
    input{
        Int flank
        File gene_sv_slop
        
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output{
        File vep_input = "vep_input.~{flank}.bed"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
    command <<<
        prep_vep_input "~{gene_sv_slop}" "vep_input.~{flank}.bed"
    >>>
}
