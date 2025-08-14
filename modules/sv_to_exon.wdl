
version 1.0

task sv_to_exon{
    input{
        Int flank

        File exon_bed
        File gene_sv_bed
        
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output{
        File exon_sv_tsv = "exon_sv.${flank}.tsv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        bedtools intersect -wo -a ~{exon_bed} -b ~{gene_sv_bed} | \
        awk '$4==$14{OFS="\t";print $4,$12,$5,$6,$7,$8,$15}' | \
        sort -k1,1 -k2,2 -k3,3n > exon_sv.~{flank}.unprocessed_info.tsv
        
        extract_SV_exon_info \
        --input exon_sv.~{flank}.unprocessed_info.tsv \
        --output exon_sv.~{flank}.tsv
    >>>
}
