process run_vep{
    tag "${sample_id}-${sex}-${hap_name}-${type}-${gene_rank}"
    publishDir "${params.intermediate_outputs_dir}", mode: 'copy'
    label "run_vep"


    input:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(type), val(gene_rank), path(converted_dipcall_vcf)

    output:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(type), val(gene_rank), path("${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}_vep115.vcf"), emit: vep115_raw_vcf

    script:
    """
    output_vcf="${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}_vep115.vcf"

    # Optional hard skip for debugging or partial pipeline runs.
    if [[ "${params.skip_vep}" == "true" ]]; then
        if [[ -s "${converted_dipcall_vcf}" ]]; then
            cp "${converted_dipcall_vcf}" "\${output_vcf}"
        else
            printf "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n" > "\${output_vcf}"
        fi
        exit 0
    fi

    # Gracefully bypass VEP if input is empty/invalid or has no variant records.
    # This keeps downstream aggregation running for low-signal haplotypes.
    if [[ "${params.skip_vep_on_empty_vcf}" == "true" ]]; then
        if [[ ! -s "${converted_dipcall_vcf}" ]] || ! grep -q '^#CHROM' "${converted_dipcall_vcf}"; then
            printf "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n" > "\${output_vcf}"
            exit 0
        fi

        if ! awk 'BEGIN{found=0} !/^#/ {found=1; exit} END{exit(found?0:1)}' "${converted_dipcall_vcf}"; then
            cp "${converted_dipcall_vcf}" "\${output_vcf}"
            exit 0
        fi
    fi

    vep -i "${converted_dipcall_vcf}" \
        --force_overwrite --vcf --buffer_size 50000 --species homo_sapiens \
        --fork ${task.cpus} -o "\${output_vcf}" --cache --offline --dir_cache ${params.cache_directory} --canonical \
        --symbol --numbers --assembly GRCh38 --use_given_ref --pick_allele_gene --pick_order biotype,rank,mane_select,mane_plus_clinical --domains --pubmed --gene_phenotype \
        --sift b --polyphen b --regulatory --total_length --af --max_af --af_1kg --custom_multi_allelic \
        --fasta ${params.OPN1MW_reference} --hgvsg --hgvs --hgvsp_use_prediction --dir_plugins ${params.plugin_dir} --plugin SpliceVault,file=${params.SPLICEVAULT} \
        --plugin Enformer,file=${params.ENFORMER} --plugin OpenTargets,file=${params.OPENTARGETS} --plugin DosageSensitivity,file=${params.DOSAGE} \
        --plugin AlphaMissense,file=${params.ALPHAMISSENSE} --plugin CADD,${params.CADD} --plugin SpliceAI,snv=${params.SPLICEAISNV},indel=${params.SPLICEAIINDEL} \
        --custom file=${params.GNOMAD},short_name=gnomADg,format=vcf,type=exact,coords=0,fields=AF \
        --custom file=${params.CLINVAR},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN

    """

}