version 1.0

task combine_annotations_ABC {
    input {
        Int flank
        Int min_support_tissue

        Float zscore_threshold

        File sv_VCF
        File genotype_VCF
        File gene_list
        File sv_gene_pairs
        File expression_file
        File maf_file
        File length_file
        File CN_file

        String annotations_dir
        String expression_field
        String expression_id_field
        String maf_mode
        String maf_field
        String length_mode
        String length_field
        String CN_mode
        String collapse_mode

        Boolean remove_control_genes
    }

    output {
        File combined_dataset = "combined_dataset.tsv"
    }

    runtime {
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        if(${remove_control_genes}){
            python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
                --vcf ${sv_VCF} \
                --genotypes ${genotype_VCF} \
                --genes ${genes_list} \
                --gene-sv ${sv_gene_pairs} \
                --annotation_dir ${annotations_dir} \
                --outfile combined_dataset.tsv \
                --expressions ${expression_file} \
                --expression-field ${expression_field} \
                --expression-id-field ${expression_id_field} \
                --maf-mode ${maf_mode} \
                --maf-file ${maf_file} \
                --maf-field ${maf_field} \
                --length-mode ${length_mode} \
                --length-file ${length_file} \
                --length-field ${length_field} \
                --CN-mode ${CN_mode} \
                --CN-file ${CN_file} \
                --collapse-mode ${collapse_mode} \
                --remove-control-genes \
                --flank ${flank} \
                --zscore-threshold ${zscore_threshold} \
                --minimum-support-tissue-count ${min_support_tissue}
        }
        if(!${remove_control_genes}){
            python3.10 scripts/executable_scripts/combine_all_annotations_ABC_polars.py \
                --vcf ${sv_VCF} \
                --genotypes ${genotype_VCF} \
                --genes ${genes_list} \
                --gene-sv ${sv_gene_pairs} \
                --annotation_dir ${annotations_dir} \
                --outfile combined_dataset.tsv \
                --expressions ${expression_file} \
                --expression-field ${expression_field}
                --expression-id-field ${expression_id_field} \
                --maf-mode ${maf_mode} \
                --maf-file ${maf_file} \
                --maf-field ${maf_field} \
                --length-mode ${length_mode} \
                --length-file ${length_file} \
                --length-field ${length_field} \
                --CN-mode ${CN_mode} \
                --CN-file ${CN_file} \
                --collapse-mode ${collapse_mode} \
                --no-remove-control-genes \
                --flank ${flank} \
                --zscore-threshold ${zscore_threshold} \
                --minimum-support-tissue-count ${min_support_tissue}

        }
    >>>
}
