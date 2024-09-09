version 1.0

task merge_enhancers {
    input {
        File enhancers
        File primary_cells
    }

    output {
        File merged_enhancers = "primary_cells_collapsed_enhancers.bed"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        python 3.10 scripts/executable_scripts/merge_enhancers.py \
            --enhancers ${enhancers} \
            --primary-cell-list ${primary_cells} \
            --out-merged-enhancers primary_cells_collapsed_enhancers.bed
    >>>
}
