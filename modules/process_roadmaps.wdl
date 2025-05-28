version 1.0

task process_roadmaps{
    input{
        File gene_sv_bed
        File roadmap_tar
        Int flank
        
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output{
        Array[File] generic_roadmap = glob("roadmap_multitissue_sv_to_gene.*.${flank}.tsv")
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
    
	String curr_roadmap_dir=basename(roadmap_tar,".tar.gz")
    
    command <<<

        tar -xzvf ~{roadmap_tar}

        split_bed_by_stateno ~{curr_roadmap_dir} "~{curr_roadmap_dir}_processed"   
        
        for i in {1..25}
        do
        sv_to_gene_roadmap ~{gene_sv_bed} "~{curr_roadmap_dir}_processed" roadmap_multitissue_sv_to_gene.${i}.~{flank}.tsv ${i}
        done

    >>>
}