version 1.0

task sv_to_gene_processing{
    input{
        Int flank
        
        File pipeline_bed
        File genes_bed

        String outdir
    }

    output{
        File gene_sv_bed = "${outdir}/gene_sv.${flank}.bed"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }


    command <<<
        bedtools intersect -a ${pipeline_bed} -b ${genes_bed} -wb | awk '{{OFS="\t";print $1,$2,$3,$4,$5,$9}}' > ${outdir}/gene_sv.${flank}.bed
    >>>
}
