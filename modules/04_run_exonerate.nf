process run_exonerate {

    tag "${sample_id}-${sex}-${params.region_name}-${hap_name}"
    publishDir "${params.intermediate_outputs_dir}", mode: 'copy'
    label "run_exonerate"


    input:
    tuple val(sample_id), val(sex), path(assembly_fa), val(hap_name)

    output:
    tuple val(sample_id), val(sex), val(hap_name), path("${sample_id}-${sex}-${params.region_name}-${hap_name}.exonerate.vulgar"), emit: exonerate_output_vulgar

    script:
    """
    exonerate --model est2genome --bestn 10 --showalignment no --showvulgar yes \
        --percent 98 \
        --ryo ">%qi|%ti|%qab-%qae|%tab-%tae|%s\\n%tas\\n" \
        --query "${params.opsin_reference}" --target "${assembly_fa}" > "${sample_id}-${sex}-${params.region_name}-${hap_name}.exonerate" 
    
    grep "^vulgar:" "${sample_id}-${sex}-${params.region_name}-${hap_name}.exonerate" > "${sample_id}-${sex}-${params.region_name}-${hap_name}.exonerate.vulgar" 2>/dev/null || true
    
    """
}