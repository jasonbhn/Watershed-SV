version 1.0

task sv_to_gene_tad{
    input{
        File TADs_dir
        File genome_bound_file
        File gene_sv_bed
    }

    output{
        File sv_to_genes_tad = "TAD_boundary_by_genes_SV.dist.${flank}.tsv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }


    command <<<
        
        bash scripts/executable_scripts/clean_TADs.sh ${TADs_dir} TAD_cleaned
        
        bedtools slop -i ${gene_sv_bed} -g ${genome_bound_file} -b 5000 > TAD_5000_flank_gene_SV.bed

        bedtools intersect -wa -wb -a TAD_5000_flank_gene_SV.bed -b TAD_cleaned/*  -filenames | \
        sort -k4,4 -k6,6 | \
        bedtools groupby -i stdin -g 4,6 -c 7 -o count_distinct | \
        awk 'BEGIN{{print "SV\tGene\tnum_TADs"}};{{OFS="\t";print}}' > TAD_boundary_by_genes_SV.dist.${flank}.tsv

    >>>
}