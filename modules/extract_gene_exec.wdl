version 1.0

task extract_gene_exec {
    input {
        File gencode_genes 
        File genome_bound_file
    }

    output {
        File genes = "genes.bed"
        File exons = "exons.bed"
        File gene_tss = "gene_tss.tsv"
        File gene_tes = "gene_tes.tsv"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        python3.10 scripts/executable_scripts/extract_gene_exec.py \
            --gencode-annotations ${gencode_genes} \
            --out-gene-bed genes.bed \
            --out-exon-bed exons.bed \
            --out-gene-tss gene_tss.tsv \
            --out-gene-tes gene_tes.tsv \
            --genome-bound ${genome_bound_file}
    >>>
}
