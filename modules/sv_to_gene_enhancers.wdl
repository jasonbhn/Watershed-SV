version 1.0

task sv_to_gene_enhancers {
    input {
        Int flank

        File merged_enhancers
        File gene_sv_bed

        String outdir
    }

    output {
        File gene_sv_enhancer = "${outdir}/enhancers_by_genes_SV.dist.${flank}.tsv"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        bedtools intersect -wa -wb -a ${merged_enhancers} -b ${gene_sv_bed} | \
        awk '{OFS="\t"; print $8,$10,$4}' | \
        sort -k1,2 -k2,2 | \
        bedtools groupby -i stdin -g 1,2 -c 3 -o max | \
        awk 'BEGIN{print "SV\tGene\tnum_enhancers_cell_types"};{OFS="\t";print $0}' \
        > ${outdir}/enhancers_by_genes_SV.dist.${flank}.tsv
    >>>
}
