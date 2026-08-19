process extract_annotation_sequences {
    tag "${sample_id}-${sex}-${hap_name}-${type}-${gene_rank}"
    publishDir "${params.intermediate_outputs_dir}", mode: 'copy'
    label "extract_annotation_sequences"

    input:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(start_position), val(end_position), val(type), val(gene_rank), path(assembly_fa)

    output:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(type), val(gene_rank), path("${sample_id}-${sex}-${hap_name}-${contig}-*-*.fa"), emit: extracted_fasta

    script:
    """
    samtools faidx "${assembly_fa}" "${contig}:${start_position}-${end_position}" > "${sample_id}-${sex}-${hap_name}-${contig}-${type}-${gene_rank}.fa"
    """
}