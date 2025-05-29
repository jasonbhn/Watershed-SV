
version 1.0

task sv_to_gene_dist {
    input {
        Int flank

        File gene_bed
        File gene_sv_slop_bed
        File gene_tss
        File gene_tes
        
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output {
        File sv_dist_to_gene = "SV_dist_to_gene.dist.${flank}.tsv"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        sv_to_gene_dist \
            --flank ~{flank} \
            --gene ~{gene_bed} \
            --gene-sv ~{gene_sv_slop_bed} \
            --in-gene-tss ~{gene_tss} \
            --in-gene-tes ~{gene_tes} \
            --out-gene-sv-dist SV_dist_to_gene.dist.~{flank}.tsv
    >>>
}
