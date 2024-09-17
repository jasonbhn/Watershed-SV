version 1.0 

task extract_rare_variants{
    input{
        File input_vcf
        File metadata

        String pipeline
        Array[String] filters
        String outdir

        Boolean filter_ethnicity
        Boolean filter_rare
    }

    output{
        File pipeline_input = "${outdir}/pipeline_input.bed"
        File vep_input = "${outdir}/vep_input.tsv"
        File pipeline_maf = "${outdir}/pipeline_maf.tsv"
        File pipeline_input_genotypes = "${outdir}/pipeline_input_genotypes.tsv"
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
            if [ "$filter_ethnicity" == "true" ]; then
              if [ "$filter_rare" == "true" ]; then 
                python3.10 scripts/executable_scripts/extract_rare_variants.py \
                --vcf $input_vcf \
                --lifted-coord $liftover_bed \
                --extract-genotype \
                --infer-rareness \
                --filter-ethnicity \
                --metadata $metadata \
                --genotype-filters '~{sep=" " filters}' \
                --out-annotsv ${outdir}/vep_input.tsv \
                --out-generic ${outdir}/pipeline_input.bed \
                --out-maf ${outdir}/pipeline_maf.tsv \
                --out-genotype ${outdir}/pipeline_input_genotypes.tsv
              else 
                python3.10 scripts/executable_scripts/extract_rare_variants.py \
                --vcf $input_vcf \
                --lifted-coord $liftover_bed \
                --extract-genotype \
                --filter-ethnicity \
                --metadata $metadata \
                --genotype-filters '~{sep=" " filters}' \
                --out-annotsv ${outdir}/vep_input.tsv \
                --out-generic ${outdir}/pipeline_input.bed \
                --out-maf ${outdir}/pipeline_maf.tsv \
                --out-genotype ${outdir}/pipeline_input_genotypes.tsv
              fi
            else
              if [ "$filter_rare" == "true" ]; then 
                python3.10 scripts/executable_scripts/extract_rare_variants.py \
                --vcf $input_vcf \
                --lifted-coord $liftover_bed \
                --extract-genotype \
                --infer-rareness \
                --genotype-filters '~{sep=" " filters}' \
                --out-annotsv ${outdir}/vep_input.tsv \
                --out-generic ${outdir}/pipeline_input.bed \
                --out-maf ${outdir}/pipeline_maf.tsv \
                --out-genotype ${outdir}/pipeline_input_genotypes.tsv
              else
                python3.10 scripts/executable_scripts/extract_rare_variants.py \
                --vcf $input_vcf \
                --lifted-coord $liftover_bed \
                --extract-genotype \
                --genotype-filters '~{sep=" " filters}' \
                --out-annotsv ${outdir}/vep_input.tsv \
                --out-generic ${outdir}/pipeline_input.bed \
                --out-maf ${outdir}/pipeline_maf.tsv \
                --out-genotype ${outdir}/pipeline_input_genotypes.tsv
              fi
            fi
          fi
        elif [ "$pipeline" == "smallset" ]; then
          if [ ! -f "pipeline_input.bed" ]; then
          python3.10 scripts/executable_scripts/extract_rare_variants.py \
          --vcf $input_vcf \
          --extract-genotype \
          --genotype-filters ${filters} \
          --out-annotsv ${outdir}/vep_input.tsv \
          --out-generic ${outdir}/pipeline_input.bed \
          --out-genotype ${outdir}/pipeline_input_genotypes.tsv
          fi     
        fi
    >>>
}
