version 1.0

task merge_enhancers {
    input {
        Int flank

        File enhancers
        File primary_cells
        File gene_sv_bed
        
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output {
        File merged_enhancers = "primary_cells_collapsed_enhancers.bed"
        File enhancers_by_genes = "enhancers_by_genes_SV.dist.${flank}.tsv"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        merge_enhancers \
            --enhancers ~{enhancers} \
            --primary-cell-list ~{primary_cells} \
            --out-merged-enhancers primary_cells_collapsed_enhancers.bed

        bedtools intersect -wa -wb -a primary_cells_collapsed_enhancers.bed -b ~{gene_sv_bed} | \
        awk '{OFS="\t";print $8,$10,$4}' | \
        sort -k1,1 -k2,2 | \
        bedtools groupby -i stdin -g 1,2 -c 3 -o max | \
        awk 'BEGIN{print "SV\tGene\tnum_enhancers_cell_types"};{OFS="\t";print $0}' \
        > enhancers_by_genes_SV.dist.~{flank}.tsv
    >>>
}