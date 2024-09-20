version 1.0

task train_test_predict{
    input{
        String training
        String testings
        String mode
        String out_prefix

        Int min_af_impute_mode
        Int min_af_value
    }

    output{
        File training_out = "${out_prefix}_training_data.tsv"
        File testing_out = "${out_prefix}_testing_data.tsv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        combine_all_annotations_ABC_polars \
        --training ${training} \
        --testings ${testings} \
        --mode ${mode} \
        --min_af_impute_mode ${min_af_impute_mode} \
        --out_prefix ${out_prefix}
    >>>
}