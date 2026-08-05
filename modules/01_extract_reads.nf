process extract_reads {

    tag "${sample_id}-${sex}-${params.region_name}"
    publishDir "${params.output_dir}", mode: 'copy'
    label "extract_reads"

    input:
    tuple val(sample_id), val(sex), path(aligned_bam), path(aligned_bam_index)

    output:
    tuple val(sample_id), val(sex), path("${sample_id}.${sex}.fastq"), emit: fastq

    script:
    """
    samtools view -@ ${task.cpus} -b ${params.samtools_flag} "${aligned_bam}" "${params.region}" > "${sample_id}.${sex}.temp.bam"


    samtools fastq "${sample_id}.${sex}.temp.bam" > "${sample_id}.${sex}.fastq"

    rm "${sample_id}.${sex}.temp.bam"
    """

}