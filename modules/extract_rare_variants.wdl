version 1.0 

task extract_rare_variants{
    input{
        File input_vcf
        File metadata

        String pipeline

        Boolean filter_ethnicity
        Boolean filter_rare
    }

    output{
        File "pipeline_input.bed"
        File "vep_input.tsv"
        File "pipeline_maf.tsv"
        File "pipeline_input_genotypes.tsv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }

    command <<<
        if [ "$pipeline" == "population" ]; then
          if [ ! -f "pipeline_input.bed" ]; then
            if [ "$filter_ethnicity" == "True" ]; then
              if [ "$filter_rare" == "True" ]; then 
                python3.10 scripts/executable_scripts/extract_rare_variants.py \
                --vcf $input_vcf \
                --lifted-coord $liftover_bed \
                --extract-genotype \
                --infer-rareness \
                --filter-ethnicity \
                --metadata $metadata \
                --genotype-filters $filters \
                --out-annotsv vep_input.tsv \
                --out-generic pipeline_input.bed \
                --out-maf pipeline_maf.tsv \
                --out-genotype pipeline_input_genotypes.tsv
              else 
                python3.10 scripts/executable_scripts/extract_rare_variants.py \
                --vcf $input_vcf \
                --lifted-coord $liftover_bed \
                --extract-genotype \
                --filter-ethnicity \
                --metadata $metadata \
                --genotype-filters $filters \
                --out-annotsv vep_input.tsv \
                --out-generic pipeline_input.bed \
                --out-maf pipeline_maf.tsv \
                --out-genotype pipeline_input_genotypes.tsv
              fi
            else
              if [ "$filter_rare" == "True" ]; then 
                python3.10 scripts/executable_scripts/extract_rare_variants.py \
                --vcf $input_vcf \
                --lifted-coord $liftover_bed \
                --extract-genotype \
                --infer-rareness \
                --genotype-filters $filters \
                --out-annotsv vep_input.tsv \
                --out-generic pipeline_input.bed \
                --out-maf pipeline_maf.tsv \
                --out-genotype pipeline_input_genotypes.tsv
              else
                python3.10 scripts/executable_scripts/extract_rare_variants.py \
                --vcf $input_vcf \
                --lifted-coord $liftover_bed \
                --extract-genotype \
                --genotype-filters $filters \
                --out-annotsv vep_input.tsv \
                --out-generic pipeline_input.bed \
                --out-maf pipeline_maf.tsv \
                --out-genotype pipeline_input_genotypes.tsv
              fi
            fi
          fi
        elif [ "$pipeline" == "smallset" ]; then
          if [ ! -f "pipeline_input.bed" ]; then
          python3.10 scripts/executable_scripts/extract_rare_variants.py \
          --vcf $input_vcf \
          --extract-genotype \
          --genotype-filters ${filters} \
          --out-annotsv vep_input.tsv \
          --out-generic pipeline_input.bed \
          --out-genotype pipeline_input_genotypes.tsv
          fi     
        fi
    >>>
}