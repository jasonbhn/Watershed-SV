version 1.0

task sv_to_geneCPG{
    input{
        File cpg_file
        File gene_sv_bed

        Int flank

        String outdir

    }

    outputs{
        File sv_to_gene_cpg_dist = "${outdir}/sv_to_gene_cpg.dist.${flank}.tsv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        awk '{{FS="\t";OFS="\t";print $2,$3,$4,$10}}' ${cpg_file} > cpgtmp.bed

        bedtools intersect -wa -wb -a ${gene_sv_bed} -b cpgtmp.bed > cpg_by_genes_SV.dist.${flank}.bed

        python3.10 scripts/executable_scripts/sv_to_gene_cpg.py --gene-sv-cpg cpg_by_genes_SV.dist.${flank}.bed --out-gene-sv-cpg ${outdir}/sv_to_gene_cpg.dist.${flank}.tsv
    >>>
}
