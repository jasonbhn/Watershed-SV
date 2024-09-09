version 1.0

task vep {
    input{
        Int flank
        String vep_cache_dir
    }

    output{
        File vep = "vep_out.csv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
}

command <<<
    python3.10 scripts/executable_scripts/prep_vep_input.py gene_sv_slop.${flank}.bed vep_input.${flank}.bed

    sort -k1,1 -k2,2n vep_input.${flank}.bed | vep \
    -o tmp.tsv \
    --format ensembl \
    --verbose \
    --cache \
    --dir ${vep_cache_dir} \
    --tab \
    --fields "Uploaded_variation,Gene,Feature_type,Consequence,IMPACT" \
    --fork 4 \
    --force \
    --regulatory \
    --overlaps \
    --distance 10000

    python3.10 scripts/executable_scripts/extract_sv_vep_annotations.py tmp.csv vep_out.csv
>>>