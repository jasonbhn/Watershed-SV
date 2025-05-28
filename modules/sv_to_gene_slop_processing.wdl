
version 1.0

task sv_to_gene_slop_processing{
    input{
        Int flank

        File genome_bound_file
        File genes_bed
        File pipeline_bed

        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output{
        File sv_gene_slop_bed = "gene_sv_slop.${flank}.bed"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        bedtools slop -g ~{genome_bound_file} -i ~{genes_bed} -b ~{flank} | \
        bedtools intersect -a ~{pipeline_bed} -b stdin -wb | \
        awk '{{OFS="\t";print $1,$2,$3,$4,$5,$9}}' > gene_sv_slop.~{flank}.bed
    >>>

}
