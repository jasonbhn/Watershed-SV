version 1.0

task process_roadmaps{
    input{
        File gene_sv_bed
        String roadmap_dir
        Int flank
    }

    output{
        Array[File] generic_roadmap = glob("roadmap_multitissue_sv_to_gene.generic.*.tsv")
        Array[File] tss_roadmap = glob("roadmap_multitissue_sv_to_gene.tss.*.tsv")
        Array[File] tes_roadmap = glob("roadmap_multitissue_sv_to_gene.tes.*.tsv")
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
        scripts/executable_scripts/sv_to_gene_roadmap.sh ${gene_sv_bed} processed_roadmaps roadmap_multitissue_sv_to_gene.generic.${i}.tsv ${i}
        scripts/executable_scripts/sv_to_gene_roadmap.sh ${gene_sv_bed} processed_roadmaps roadmap_multitissue_sv_to_gene.tss.${i}.tsv ${i}
        scripts/executable_scripts/sv_to_gene_roadmap.sh ${gene_sv_bed} processed_roadmaps roadmap_multitissue_sv_to_gene.tes.${i}.tsv ${i}
        done

    >>>
}
