
version 1.0

task sv_to_gene_tes_flank_processing {
    input {
        Int flank

        File genome_bound_file
        File genes_bed
        File pipeline_bed

        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output {
        File gene_sv_flank_bed = "gene_sv.${flank}.bed" 
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        bedtools flank -g ~{genome_bound_file} -i ~{genes_bed} -r ~{flank} -l 0 -s | \
        bedtools intersect -a ~{pipeline_bed} -b stdin -wb | \
        awk '{{OFS="\t";print $1,$2,$3,$4,$5,$9}}' > gene_sv.~{flank}.bed
    >>>
}
