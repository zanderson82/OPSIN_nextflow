process convert_gff_to_bed {

    tag "${sample_id}-${sex}-${params.region_name}-${hap_name}"
    label "convert_gff_to_bed"

    input:
    tuple val(sample_id), val(sex), val(hap_name), path(vulgar_to_gff_output)

    output:
    tuple val(sample_id), val(sex), val(hap_name), path("${sample_id}_${hap_name}.bed"), emit: gff_to_bed_output

    script:
    """
    gff_to_bed2.py \
        --gff_file "${vulgar_to_gff_output}" \
        --sample_name "${sample_id}" \
        --sex "${sex}" \
        --hap_name "${hap_name}" \
        --output_dir "."
    """
}