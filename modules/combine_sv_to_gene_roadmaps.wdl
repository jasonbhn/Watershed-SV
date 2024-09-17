version 1.0

task combine_sv_to_gene_roadmaps {
    input {
        Int flank

        String roadmap_dir
        String outdir
    }

    output {
        File gene_sv_roadmaps = "${outdir}/combined_roadmaps.dist.${flank}.tsv"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        python3.10 scripts/executable_scripts/combine_roadmaps.py \
            --gene-sv-roadmap-dir ${roadmap_dir} \
            --out-combined-roadmap ${outdir}/combined_roadmaps.dist.${flank}.tsv
    >>>
}
