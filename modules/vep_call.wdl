
version 1.0

task vep_call {
    input{
        File vep_cache_tar
        File vep_input
        String docker
        Int memory
        Int disk_space
        Int ncpu
    }

    output{
        File vep_out = "tmp.tsv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
    
    command <<<
        mkdir -p "vep_cache"
        tar -xzvf ~{vep_cache_tar} -C "vep_cache/"
        ls
        ls vep_cache
        sort -k1,1 -k2,2n "~{vep_input}" | vep \
        -o "tmp.tsv" \
        --format ensembl \
        --verbose \
        --cache \
        --offline \
        --dir_cache vep_cache/vep_cache \
        --tab \
        --fields "Uploaded_variation,Gene,Feature_type,Consequence,IMPACT" \
        --fork 4 \
        --force \
        --regulatory \
        --overlaps \
        --distance 10000
    >>>
}
