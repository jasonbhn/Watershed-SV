input{
    File gene_sv_bed
    File ABC_enhancers

    String output_prefix
}

output{
    File sv_to_gene_flank = "${output_prefix}.${flank}.tsv"
}

runtime{
    docker: "${docker}"
    memory: "${memory}GB"
    disks: "local-disk ${disk_space} HDD"
    cpu: "${ncpu}"
}

command <<<
    
    bedtools intersect -wa -wb -a ${gene_sv_bed} -b ${ABC_enhancers} | \
    awk 'BEGIN{print "SV\tGene\tis_ABC_SV"}; split($6,a,".") {OFS="\t";if (a[1]==$17) print $4,a[1],1}' > ${output_prefix}.${flank}.tsv

>>>