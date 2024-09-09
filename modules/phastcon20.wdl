version 1.0

task phastcon20{
    input{
        File gene_sv_bed
        File phcon20_file

        String outut_prefix
    }

    output{
        File ${output_prefix}.${flank}.tsv
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
    --in-bigwig ${phcon20_file} \
    --bigwig-name "top10_phastCON" \
    --stat-method "top10_mean" \
    --score-upper-limit 1 \
    --score-lower-limit 0 \
    --out-gene-sv-score ${output_prefix}.${flank}.tsv
>>>