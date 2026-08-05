process run_hifiasm_XX {

    tag "${sample_id}-${sex}-${params.region_name}"
    publishDir "${params.output_dir}", mode: 'copy'
    label "run_hifiasm_XX"

    input:
    tuple val(sample_id), val(sex), path(fastq), val(hg_size)

    output:
    tuple val(sample_id), val(sex), path("${sample_id}.${sex}.${params.region_name}.asm.hap1.p_ctg.fa"), optional: true, emit: hap1_fa
    tuple val(sample_id), val(sex), path("${sample_id}.${sex}.${params.region_name}.asm.hap2.p_ctg.fa"), optional: true, emit: hap2_fa
    tuple val(sample_id), val(sex), path("${sample_id}.${sex}.${params.region_name}.asm.hap1.p_ctg.fa.fai"), optional: true, emit: hap1_fai
    tuple val(sample_id), val(sex), path("${sample_id}.${sex}.${params.region_name}.asm.hap2.p_ctg.fa.fai"), optional: true, emit: hap2_fai

    script:
    """
    ${projectDir}/resources/hifiasm/hifiasm -o "${sample_id}.${sex}.${params.region_name}.asm" --ont -t${task.cpus} --hg-size "${hg_size}" "${fastq}"

    awk '/^S/{print ">"\$2"\\n"\$3}' "${sample_id}.${sex}.${params.region_name}.asm.bp.hap1.p_ctg.gfa" > "${sample_id}.${sex}.${params.region_name}.asm.hap1.p_ctg.fa"
    if [ -s "${sample_id}.${sex}.${params.region_name}.asm.hap1.p_ctg.fa" ]; then
        samtools faidx "${sample_id}.${sex}.${params.region_name}.asm.hap1.p_ctg.fa"
    else
        rm -f "${sample_id}.${sex}.${params.region_name}.asm.hap1.p_ctg.fa"
    fi

    awk '/^S/{print ">"\$2"\\n"\$3}' "${sample_id}.${sex}.${params.region_name}.asm.bp.hap2.p_ctg.gfa" > "${sample_id}.${sex}.${params.region_name}.asm.hap2.p_ctg.fa"
    if [ -s "${sample_id}.${sex}.${params.region_name}.asm.hap2.p_ctg.fa" ]; then
        samtools faidx "${sample_id}.${sex}.${params.region_name}.asm.hap2.p_ctg.fa"
    else
        rm -f "${sample_id}.${sex}.${params.region_name}.asm.hap2.p_ctg.fa"
    fi

    awk '/^S/{print ">"\$2"\\n"\$3}' "${sample_id}.${sex}.${params.region_name}.asm.bp.p_ctg.gfa" > "${sample_id}.${sex}.${params.region_name}.asm.p_ctg.fa"
    samtools faidx "${sample_id}.${sex}.${params.region_name}.asm.p_ctg.fa"
    """

}