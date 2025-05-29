
version 1.0

task TAD_annotation_clean_up{
    input{
        File TADs_tar
        
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output{
        Array[File] cleaned_TAD_beds = glob("*.bed")
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    String curr_tad_dir = basename(TADs_tar,".tar.gz")
    command <<<
        tar -xzvf ~{TADs_tar}
        clean_TADs_WDL ~{curr_tad_dir}
    >>>
}
