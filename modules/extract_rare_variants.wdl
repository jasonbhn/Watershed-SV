version 1.0 

task extract_rare_variants{
    input{
        File input_vcf
        File metadata

        String pipeline
        Array[String] filters
        

        Boolean filter_ethnicity
        Boolean filter_rare
        
        String docker
        Int memory
        Int disk_space
        Int ncpu
        
        File? liftover_bed

    }
    
	output{
        File pipeline_input = "pipeline_input.bed"
        File pipeline_input_genotypes = "pipeline_input_genotypes.tsv"
    }

    runtime{
        docker: "${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"
    }
    command <<<
        if [ "~{pipeline}" == "population" ]; then
          if [ ! -f "pipeline_input.bed" ]; then
            if [ "~{filter_ethnicity}" == "true" ]; then
              if [ "~{filter_rare}" == "true" ]; then
                extract_rare_variants \
                --vcf "~{input_vcf}" \
                ~{"--lifted-coord " + liftover_bed} \
                --extract-genotype \
                --infer-rareness \
                --filter-rare \
                --filter-ethnicity \
                --metadata "~{metadata}" \
                --genotype-filters ~{sep=" " filters} \
                --out-annotsv vep_input.tsv \
                --out-generic pipeline_input.bed \
                --out-maf pipeline_maf.tsv \
                --out-genotype pipeline_input_genotypes.tsv
              else 
                extract_rare_variants \
                --vcf ~{input_vcf} \
                ~{"--lifted-coord " + liftover_bed} \
                --infer-rareness \
                --no-filter-rare \
                --extract-genotype \
                --filter-ethnicity \
                --metadata ~{metadata} \
                --genotype-filters ~{sep=" " filters} \
                --out-annotsv vep_input.tsv \
                --out-generic pipeline_input.bed \
                --out-maf pipeline_maf.tsv \
                --out-genotype pipeline_input_genotypes.tsv
              fi
            else
              if [ "~{filter_rare}" == "true" ]; then 
                extract_rare_variants \
                --vcf ~{input_vcf} \
                ~{"--lifted-coord " + liftover_bed} \
                --extract-genotype \
                --no-filter-ethnicity \
                --infer-rareness \
                --filter-rare \
                --genotype-filters ~{sep=" " filters} \
                --out-annotsv vep_input.tsv \
                --out-generic pipeline_input.bed \
                --out-maf pipeline_maf.tsv \
                --out-genotype pipeline_input_genotypes.tsv
              else
                extract_rare_variants \
                --vcf ~{input_vcf} \
                ~{"--lifted-coord " + liftover_bed} \
                --no-filter-ethnicity \
                --infer-rareness \
                --no-filter-rare \
                --extract-genotype \
                --genotype-filters ~{sep=" " filters} \
                --out-annotsv vep_input.tsv \
                --out-generic pipeline_input.bed \
                --out-maf pipeline_maf.tsv \
                --out-genotype pipeline_input_genotypes.tsv
              fi
            fi
          fi
        elif [ "~{pipeline}" == "smallset" ]; then
          if [ ! -f "pipeline_input.bed" ]; then
          extract_rare_variants \
          --vcf ~{input_vcf} \
          --extract-genotype \
          --no-infer-rareness \
          --no-filter-ethnicity \
          --no-filter-rare \
          --genotype-filters ~{sep=" " filters} \
          --out-annotsv vep_input.tsv \
          --out-generic pipeline_input.bed \
          --out-genotype pipeline_input_genotypes.tsv
          fi     
        fi
    >>>
}