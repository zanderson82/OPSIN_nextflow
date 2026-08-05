process run_hifiasm_XY {

    tag "${sample_id}-${sex}-${params.region_name}"
    publishDir "${params.output_dir}", mode: 'copy'
    label "run_hifiasm_XY"

    input:
    tuple val(sample_id), val(sex), path(fastq), val(hg_size)
    

    output:
    tuple val(sample_id), val(sex), path("${sample_id}.${sex}.${params.region_name}.asm.p_ctg.fa"), emit: primary_fa

    script:
    """
    hifiasm -o "${sample_id}.${sex}.${params.region_name}.asm" --ont -t${task.cpus} ${params.hifiasm_opts_XY} --hg-size "${hg_size}" "${fastq}"

    awk '/^S/{print ">"\$2"\\n"\$3}' "${sample_id}.${sex}.${params.region_name}.asm.bp.p_ctg.gfa" > "${sample_id}.${sex}.${params.region_name}.asm.p_ctg.fa"
    samtools faidx "${sample_id}.${sex}.${params.region_name}.asm.p_ctg.fa"
    """

}