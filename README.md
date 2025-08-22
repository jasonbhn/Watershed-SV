# WDL input description

* `File input_vcf`
  * SV VCF file. Must contain INFO/MODECN and FORMAT/CN fields for CNVs, and a INFO/SVLEN field
* `File metadata`: Tab-delimited file containing SubjectID (must match some of the IDs of `input_vcf`), and a `RACE` column used if subsetting by population group.
  * Only used if you set `filter_ethnicity`=True, but you are still required to provide a file.
* `String pipeline`: "population" or "smallset". “If data is of sufficient size, ie > 100, select population.” "smallset" disables `filter_ethnicity` and `filter_rare` options.
* `Array[String] filters`: Keeps only the variants which has one of these strings in the input_vcf's FILTER column. By default, ["PASS", “.”, “MATCH_1KG”].
* `Boolean filter_ethnicity`: if True, filter for EUR ancestry.
* `Boolean filter_rare`: if True, filter for variants w/ MAF < 0.01.
* `File gencode_genes`: GENCODE gtf file
* `File genome_bound_file`: Tab-delimited file with 1st col chr/contig name and 2nd col the length in basepairs. Analysis will be restricted to the chrs/contigs listed.
* `Int flank`: “how much flanking up and downstream of genes to consider. Usually use 100000 or 10000.”
* `File crm_remap`: REMAP CRM annotation file
* `File ABC_enhancers`: Activity-By-Contact annotation file
* `File bw_GC`: GC content annotation file
* `File bw_linsight`: LINSIGHT scores annotation file
* `File bw_CADD`: CADD scores annotation file
* `File bw_phastcon`: PhastCons20way scores annotation file
* `File cpg_file`: CpG islands annotation file
* `File TADs_tar`: Archive of TAD annotations
* `File enhancers`: Enhancers annotation file
* `File primary_cells`: File listing of primary cell types
* `File roadmap_tar`: Archive of ROADMAP annotations
* `File vep_cache_tar`: VEP cache archive (homo_sapiens_merged_vep_112_GRCh38.tar.gz)
* `File gene_list`: Tab-delimited file with chrom, start, end, gene, gene_type, and Strand columns.
* `File expression_file`: Tab-delimited file of expression z-scores for each gene (3 columns: gene, Individual, and z-score)
* `String maf_mode`: "upload" a custom SV→MAF tsv, or “extract” from the input_vcf file.
* `File? maf_file`: If you set maf_mode == "upload", provide here the file name of a tab-delimited file w/ 2 columns: SV id and AF.
* `String maf_field`: If you set maf_mode == "extract", provide the field name of allele frequency from the input_vcf's INFO column (usually "AF").
* `Int min_support_tissue`: Minimum number of samples with non-NaN measurements for a tissue. If a tissue type has fewer samples than this, it will be dropped.
* `Float zscore_threshold`: Z-score threshold for outlier calling. Usually 3.
* `String expression_field`: Prefix that all the expression data colnames of the expression_file input have, like “ENSG”.
* `String expression_id_field`: Name of the subject ID column of `expression_file`.
* `String length_mode`: “upload-SV” or "upload-VNTR" to get length info from length_file input, or "extract" length info from input_vcf.
* `File? length_file`: If you set length_mode == "upload-SV", file name of tab-delimited file w/ 2 columns: SV id and SV length.
* `String length_field`: If you set length_mode == "extract", field name of AF in the input_vcf's INFO column.
* `String CN_mode`: "upload" to get copy number from CN_file, or "extract" from other information (it will be calculated based on the SVTYPE of each SV)
* `File? CN_file`: If you set CN_mode == “upload”, a Tab-delimited file with SubjectID, SV, and CN columns.
* `String collapse_mode`: "gene" to collapse all rare SVs to their respective nearby gene, or "gene-sv" to evaluate per-gene per-SV.
* `Boolean remove_control_genes`: Whether to remove genes with no outliers or not.
