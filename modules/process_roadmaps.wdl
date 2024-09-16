version 1.0

task process_roadmaps{
    input{
        String roadmap_dir
        Int flank
    }

    output{
        File generic_roadmap = "roadmap_multitissue_sv_to_gene.generic.${i}.tsv"
        File tss_roadmap = "roadmap_multitissue_sv_to_gene.tss.${i}.tsv"
        File tes_roadmap = "roadmap_multitissue_sv_to_gene.tes.${i}.tsv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<

        scripts/executable_scripts/split_bed_by_stateno.sh ${roadmap_dir} ./processed_roadmaps    
        
        for i in {1..25}
        do
        scripts/executable_scripts/sv_to_gene_roadmap.sh ${outdir}/intermediates/gene_sv.${flank}.bed processed_roadmaps roadmap_multitissue_sv_to_gene.generic.${i}.tsv ${i}
        scripts/executable_scripts/sv_to_gene_roadmap.sh gene_sv.${flank}.bed processed_roadmaps roadmap_multitissue_sv_to_gene.tss.${i}.tsv ${i}
        scripts/executable_scripts/sv_to_gene_roadmap.sh gene_sv.${flank}.bed processed_roadmaps roadmap_multitissue_sv_to_gene.tes.${i}.tsv ${i}
        done

    >>>
}