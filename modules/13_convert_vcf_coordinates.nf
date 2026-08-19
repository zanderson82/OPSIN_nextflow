process convert_vcf_coordinates {
    tag "${sample_id}-${sex}-${hap_name}-${type}-${gene_rank}"
    publishDir "${params.intermediate_outputs_dir}", mode: 'copy'
    label "convert_vcf_coordinates"


    input:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(type), val(gene_rank), path(dipcall_pair_vcf)


    output:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(type), val(gene_rank), path("${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}_converted.vcf"), emit: converted_dipcall_vcf

    script:
    """
    convert-opsin-vcf-coordinates.py "${sample_id}" "${hap_name}" "${sex}" "${gene_rank}" "${type}" "${contig}" "${dipcall_pair_vcf}" "${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}_converted.vcf"
    """

}