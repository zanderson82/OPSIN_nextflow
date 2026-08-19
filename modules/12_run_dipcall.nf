process run_dipcall {
    tag "${sample_id}-${sex}-${hap_name}-${type}-${gene_rank}"
    publishDir "${params.intermediate_outputs_dir}", mode: 'copy'
    label "run_dipcall"

    input:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(type), val(gene_rank), path(gene_reference), path(extracted_fasta)

    output:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(type), val(gene_rank), path("${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}.pair.vcf"), emit: dipcall_pair_vcf


    script:
    """
    samtools faidx "${gene_reference}"
    ${projectDir}/bin/dipcall.kit/run-dipcall -t ${task.cpus} -m "${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}" \
        "${gene_reference}" "${extracted_fasta}" "${extracted_fasta}" > "${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}.mak"

        make -j8 -f "${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}.mak"

        if [ -f "${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}.pair.vcf.gz" ]; then
            gunzip -f "${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}.pair.vcf.gz"
        fi

    """
}

