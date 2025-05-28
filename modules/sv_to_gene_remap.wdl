
version 1.0

task sv_to_gene_remap {
    input {
        Int flank

        File remap_crm
        File gene_sv_bed
                
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output {
        File remap_crm_sv_tsv = "remap_crm_sv.dist.${flank}.tsv"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
    
    command <<<
        zcat ~{remap_crm} | awk '{OFS="\t";print $1,$2,$3,$5}' | \
        bedtools intersect -wa -wb -a stdin -b ~{gene_sv_bed} | \
        awk '{OFS="\t";print $8,$10,$4}' | \
        sort -k1,1 -k2,2 | \
        bedtools groupby -i stdin -g 1,2 -c 3 -o max | \
        awk 'BEGIN{print "SV\tGene\tremap_crm_score"};{OFS="\t";print $0}' \
        > remap_crm_sv.dist.~{flank}.tsv
    >>>
}
