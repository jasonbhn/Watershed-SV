version 1.0

task combine_sv_to_gene_roadmaps {
    input {
        Int flank

        String roadmap_dir
        
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output {
        File gene_sv_roadmaps = "combined_roadmaps.dist.${flank}.tsv"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        combine_roadmaps \
            --gene-sv-roadmap-dir ~{roadmap_dir} \
            --out-combined-roadmap combined_roadmaps.dist.~{flank}.tsv
    >>>
}