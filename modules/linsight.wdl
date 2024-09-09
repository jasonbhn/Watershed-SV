version 1.0

task linsight{
    input{
        File linsight_file
        File gene_sv_bed
        
        String output_prefix
    }

    output{
        File LINSIGHT_by_genes_SV.dist.${flank}.tsv
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
}


command <<<
    python3.10 scripts/executable_scripts/sv_to_gene_bw_scores.py \
    --gene-sv ${gene_sv_bed} \
    --in-bigwig ${linsight_file} \
    --bigwig-name "top10_LINSIGHT" \
    --stat-method "top10_mean" \
    --score-upper-limit 1 \
    --score-lower-limit 0 \
    --out-gene-sv-score ${output_prefix}.${flank}.tsv
>>>