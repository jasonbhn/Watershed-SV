version 1.0

task extract_gene_exec {
    input {
        File gencode_genes 
        File genome_bound_file

        String outdir
    }

    output {
        File genes = "${outdir}/genes.bed"
        File exons = "${outdir}/exons.bed"
        File gene_tss = "${outdir}/gene_tss.tsv"
        File gene_tes = "${outdir}/gene_tes.tsv"
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
            --out-gene-bed ${outdir}/genes.bed \
            --out-exon-bed ${outdir}/exons.bed \
            --out-gene-tss ${outdir}/gene_tss.tsv \
            --out-gene-tes ${outdir}/gene_tes.tsv \
            --genome-bound ${genome_bound_file}
    >>>
}
