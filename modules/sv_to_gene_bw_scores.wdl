version 1.0

task sv_to_gene_bw_scores {
    input{
        File bw 
        File gene_sv_bed
    
        String name
        String stat_method
        String outdir
    
        Int upper_limit
        Int lower_limit
        Int flank
    }
    
    output{
        File gene_sv_score = "${outdir}/${name}_by_genes_SV.${flank}.tsv"
    }
    
    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
    
    command <<<
        python3.10 scripts/executable_scripts/sv_to_gene_bw_scores.py \
        --gene-sv ${gene_sv_bed} \
        --in-bigwig ${bw} \
        --bigwig-name ${name} \
        --stat-method ${stat_method} \
        --score-upper-limit ${upper_limit} \
        --score-lower-limit ${lower_limit} \
        --out-gene-sv-score ${outdir}/${name}_by_genes_SV.{flank}.tsv
        >>>
}
