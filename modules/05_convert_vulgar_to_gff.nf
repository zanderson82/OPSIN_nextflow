process convert_vulgar_to_gff {

    tag "${sample_id}-${sex}-${params.region_name}-${hap_name}"
    publishDir "${params.intermediate_outputs_dir}", mode: 'copy'
    label "convert_vulgar_to_gff"

    input:
    tuple val(sample_id), val(sex), val(hap_name), path(exonerate_output_vulgar)

    output:
    tuple val(sample_id), val(sex), val(hap_name), path("${sample_id}-${sex}-${params.region_name}-${hap_name}.gff"), emit: vulgar_to_gff_output

    script:
    """
    awk '{
        query_name = \$2
        target_name = \$6
        target_start = \$7
        target_end = \$8
        target_strand = \$9
        score = \$10
        
        query_lower = tolower(query_name)
        
        if (query_lower ~ /lcr/) {
            type = "LCR"
        } else if (query_lower ~ /lw/ && query_lower ~ /exon/) {
            type = "OPN1LW_exon5"
        } else if (query_lower ~ /mw/ && query_lower ~ /exon/) {
            type = "OPN1MW_exon5"
        } else {
            next
        }
        
        if (target_start > target_end) {
            temp = target_start
            target_start = target_end
            target_end = temp
        }
        
        print target_name "\\t" "exonerate:est2genome" "\\t" "similarity" "\\t" target_start "\\t" target_end "\\t" score "\\t" target_strand "\\t" "." "\\t" "query=" query_name ";score=" score
    }' "${exonerate_output_vulgar}" > "${sample_id}-${sex}-${params.region_name}-${hap_name}.gff"
    """
}