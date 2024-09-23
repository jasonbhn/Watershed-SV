version 1.0

task collect_files {
    input {
        Array[File] files
    }

    output {
        String all_files = "all_files.tar"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        tar -cvz --file=all_files.tar ~{" " $files}
    >>>
}
