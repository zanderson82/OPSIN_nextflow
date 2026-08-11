process align_to_assembly {

    tag "${sample_id}-${sex}-${params.region_name}"
    publishDir "${params.intermediate_outputs_dir}", mode: 'copy'
    label "align_to_assembly"

    input:
    tuple val(sample_id), val(sex), path(assembly_fa), path(fastq)

    output:
    tuple val(sample_id), val(sex), path("${sample_id}.${sex}.${params.region_name}.reads_to_assembly.bam"), emit: reads_to_assembly_bam
    tuple val(sample_id), val(sex), path("${sample_id}.${sex}.${params.region_name}.reads_to_assembly.bam.bai"), emit: reads_to_assembly_bai

    script:
    """
    minimap2 -t ${task.cpus} -ax map-ont -secondary=no "${assembly_fa}" "${fastq}" | \
        samtools sort -o "${sample_id}.${sex}.${params.region_name}.reads_to_assembly.bam" -
    samtools index "${sample_id}.${sex}.${params.region_name}.reads_to_assembly.bam"
    """

}