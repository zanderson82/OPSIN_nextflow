process combine_SNP_analysis {
    publishDir "${params.output_dir}", mode: 'copy'
    label "concatenate_SNP_results"

    input:
    path(tsvs)


    output:
    path("${params.final_output_name}_combined_SNP_analysis.tsv"), emit: combined_snp_tsv
    path("*_combined_SNP_analysis.tsv"), emit: per_sample_combined_tsvs

    script:
    """
    combine_SNP_analysis2.py \
    "${params.final_output_name}_combined_SNP_analysis.tsv" \
    ${tsvs}
    """

}