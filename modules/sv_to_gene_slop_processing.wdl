version 1.0

task sv_to_gene_slop_processing{
    input{
        Int flank

        File genome_bound_file
        File genes_bed
        File pipeline_bed
    }

    output{
        File "gene_sv_slop.${flank}.bed"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        if [ ! -f "gene_sv_slop.${flank}.bed" ]; then 
        bedtools slop -g ${genome_bound_file} -i ${genes_bed} -b ${flank} | 
        bedtools intersect -a ${pipeline_bed} -b stdin -wb | 
        awk '{{OFS="\t";print $1,$2,$3,$4,$5,$9}}' > gene_sv_slop.${flank}.bed
    >>>

}