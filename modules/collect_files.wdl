version 1.0

task collect_files {
    input {
        Array[File] files
        String outdir
    }

    output {
        String all_files = "$outdir"                
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        mkdir -p $outdir
        mv ~{" " $files} $outdir
    >>>
}
