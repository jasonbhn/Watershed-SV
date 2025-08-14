version 1.0

task combine_annotations_ABC {
    input {
        Int flank
        Int min_support_tissue
        Float zscore_threshold
        File collapse_instruction
        File sv_VCF
        File genotype_VCF
        File gene_list
        File expression_file
        File? maf_file
        File? length_file
        File? CN_file
        Array [File] gene_body_files
        Array [File] tss_files
        Array [File] tes_files
        String expression_field
        String expression_id_field
        String maf_mode
        String maf_field
        String length_mode
        String length_field
        String CN_mode
        String collapse_mode
        Boolean remove_control_genes
        String docker
        Int memory
        Int disk_space
        Int ncpu
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
        mkdir annotations_dir
        mkdir annotations_dir/intermediates
        mkdir annotations_dir/intermediates_tss_flank
        mkdir annotations_dir/intermediates_tes_flank
        mv -t annotations_dir/intermediates ~{sep=' ' gene_body_files}
        mv -t annotations_dir/intermediates_tss_flank ~{sep=' ' tss_files}
        mv -t annotations_dir/intermediates_tes_flank ~{sep=' ' tes_files}
        
        ls annotations_dir/intermediates
        
        if [ "~{remove_control_genes}" == "true" ]; then
            combine_all_annotations_ABC_polars \
                --collapse-annotation-instruction "~{collapse_instruction}" \
                --vcf "~{sv_VCF}" \
                --genotypes "~{genotype_VCF}" \
                --genes "~{gene_list}" \
                --gene-sv "annotations_dir/intermediates/gene_sv_slop.~{flank}.bed" \
                --annotation-dir annotations_dir \
                --outfile combined_dataset.tsv \
                --expressions "~{expression_file}" \
                --expression-field "~{expression_field}" \
                --expression-id-field "~{expression_id_field}" \
                --maf-mode "~{maf_mode}" \
                --maf-file "~{maf_file}" \
                --maf-field "~{maf_field}" \
                --length-mode "~{length_mode}" \
                --length-file "~{length_file}" \
                --length-field "~{length_field}" \
                --CN-mode "~{CN_mode}" \
                --CN-file "~{CN_file}" \
                --collapse-mode "~{collapse_mode}" \
                --remove-control-genes \
                --flank "~{flank}" \
                --zscore-threshold "~{zscore_threshold}" \
                --minimum-support-tissue-count "~{min_support_tissue}" \
                --filter-rare
        elif [ "~{remove_control_genes}" == "false" ]; then
            combine_all_annotations_ABC_polars \
                --collapse-annotation-instruction "~{collapse_instruction}" \
                --vcf "~{sv_VCF}" \
                --genotypes "~{genotype_VCF}" \
                --genes "~{gene_list}" \
                --gene-sv "annotations_dir/intermediates/gene_sv_slop.~{flank}.bed" \
                --annotation-dir annotations_dir \
                --outfile combined_dataset.tsv \
                --expressions "~{expression_file}" \
                --expression-field "~{expression_field}" \
                --expression-id-field "~{expression_id_field}" \
                --maf-mode "~{maf_mode}" \
                --maf-file "~{maf_file}" \
                --maf-field "~{maf_field}" \
                --length-mode "~{length_mode}" \
                --length-file "~{length_file}" \
                --length-field "~{length_field}" \
                --CN-mode "~{CN_mode}" \
                --CN-file "~{CN_file}" \
                --collapse-mode "~{collapse_mode}" \
                --no-remove-control-genes \
                --flank "~{flank}" \
                --zscore-threshold "~{zscore_threshold}" \
                --minimum-support-tissue-count "~{min_support_tissue}" \
                --filter-rare
        else
            echo "Invalid parameter: remove_control_gene"
        fi
    >>>
}
