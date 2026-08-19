process parse_gff {
    tag "${sample_id}-${sex}-${params.region_name}-${hap_name}"
    publishDir "${params.intermediate_outputs_dir}", mode: 'copy'
    label "parse_gff"

    input:
    tuple val(sample_id), val(sex), val(hap_name), path(vulgar_to_gff_output)

    output:
    tuple val(sample_id), val(sex), val(hap_name), path("${sample_id}-${sex}-${params.region_name}-${hap_name}_gene_coordinates.tsv"), emit: gene_coordinate_tsv



    script: 
    """
    opsin_array_gff_parse_strand_aware.py "${vulgar_to_gff_output}" "${sample_id}" "${hap_name}" "${sex}" "${sample_id}-${sex}-${params.region_name}-${hap_name}_gene_coordinates.tsv"
    """

}   