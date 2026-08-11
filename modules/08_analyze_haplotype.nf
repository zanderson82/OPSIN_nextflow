process analyze_haplotype {

    tag "${sample_id}-${sex}-${params.region_name}-${hap_name}"
    publishDir "${params.output_dir}", mode: 'copy'
    label "analyze_haplotype"

    input:
    tuple val(sample_id), val(sex), val(hap_name), path(stats_gff)

    output:
    tuple val(sample_id), val(sex), val(hap_name), path("${sample_id}-${sex}-${params.region_name}-${hap_name}.analysis.tsv"), emit: haplotype_analysis

    script:
    """
    python_opsin_processing_V6.py \
        "${stats_gff}" \
        "${sample_id}" \
        "${sex}" \
        "${hap_name}" > "${sample_id}-${sex}-${params.region_name}-${hap_name}.analysis.tsv"
    """
}
