process cat_assemblies {

    tag "${sample_id}-${sex}-${params.region_name}"
    publishDir "${params.output_dir}", mode: 'copy'
    label "cat_assemblies"

    input:
    tuple val(sample_id), val(sex), path(hap1_fa), path(hap2_fa)


    output:
    tuple val(sample_id), val(sex), path("${sample_id}.${sex}.combined-diploid.fa"), emit: combined_diploid_fa

    script:
    """
    cat "${hap1_fa}" "${hap2_fa}" > "${sample_id}.${sex}.combined-diploid.fa"
    samtools faidx "${sample_id}.${sex}.combined-diploid.fa"
    """

}