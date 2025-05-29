version 1.0

task sv_to_gene_bw_scores {
    input{
        File bw 
        File gene_sv_bed
    
        String name
        String stat_method
    
        Float upper_limit
        Float lower_limit
        Int flank
        
        String docker
        Int memory
        Int disk_space
        Int ncpu
        
    }
    
    output{
        File gene_sv_score = "${name}_by_genes_SV.${flank}.tsv"
    }
    
    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
    
    command <<<
        sv_to_gene_bw_scores \
        --gene-sv ~{gene_sv_bed} \
        --in-bigwig ~{bw} \
        --bigwig-name ~{name} \
        --stat-method ~{stat_method} \
        --score-upper-limit ~{upper_limit} \
        --score-lower-limit ~{lower_limit} \
        --out-gene-sv-score ~{name}_by_genes_SV.~{flank}.tsv
        >>>
}