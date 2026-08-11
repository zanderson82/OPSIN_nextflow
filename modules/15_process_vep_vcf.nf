process process_vep_vcf {
    tag "${sample_id}-${sex}-${hap_name}-${type}-${contig}-${gene_rank}"
    publishDir "${params.intermediate_outputs_dir}", mode: 'copy'
    label "process_vep_vcf"


    input:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(type), val(gene_rank), val(gene_ref), path(vep114_raw_vcf)

    output:
    tuple val(sample_id), val(sex), val(hap_name), val(contig), val(type), val(gene_rank), val(gene_ref), path("${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}-${gene_ref}.vep_SNP_analysis.tsv"), emit: processed_vep_tsv


    script:

    """
    opsin_vep_snp_analysis_target_codons.py \
    "${vep114_raw_vcf}" \
    "${gene_ref}" \
    "${sample_id}-${sex}-${hap_name}-${gene_rank}-${type}-${contig}-${gene_ref}.vep_SNP_analysis.tsv"
    
    """
}