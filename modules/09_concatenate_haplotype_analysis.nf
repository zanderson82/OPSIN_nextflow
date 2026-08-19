process concatenate_results {
    publishDir "${params.output_dir}", mode: 'copy'
    label "concatenate_results"

    input:
    path(tsvs)

    output:
    path("${params.final_output_name}_summary.tsv"), emit: summary_tsv

   script:
    """
    python_opsin_processing_V6.py --header > "${params.final_output_name}_summary.tsv"
    cat ${tsvs} | sort -t\$'\\t' -k1,1 >> "${params.final_output_name}_summary.tsv"
    """
}