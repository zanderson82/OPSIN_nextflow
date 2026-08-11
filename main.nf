nextflow.enable.dsl=2

// Pull modules
include { extract_reads } from "./modules/01_extract_reads.nf"
include { run_hifiasm_XY } from "./modules/02.1_run_hifiasm_XY.nf"
include { run_hifiasm_XX } from "./modules/02.2_run_hifiasm_XX.nf"
include { cat_assemblies } from "./modules/02.3_cat_assemblies.nf"
include { align_to_assembly } from "./modules/03_align_to_assembly.nf"
include { run_exonerate } from "./modules/04_run_exonerate.nf"
include { convert_vulgar_to_gff } from "./modules/05_convert_vulgar_to_gff.nf"
include { convert_gff_to_bed } from "./modules/06_convert_gff_to_bed.nf"
include { add_alignment_stats } from "./modules/07_add_alignment_stats.nf"
include { analyze_haplotype } from "./modules/08_analyze_haplotype.nf"
include { concatenate_results } from "./modules/09_concatenate_haplotype_analysis.nf"
include { parse_gff } from "./modules/10_parse_gff_annotation_coordinates.nf"
include { extract_annotation_sequences } from "./modules/11_extract_annotation_sequences.nf"
include { run_dipcall } from "./modules/12_run_dipcall.nf"
include { convert_vcf_coordinates } from "./modules/13_convert_vcf_coordinates.nf"
include { run_vep } from "./modules/14_run_vep.nf"
include { process_vep_vcf } from "./modules/15_process_vep_vcf.nf"
include { combine_SNP_analysis } from "./modules/16_combine_SNP_analysis_outputs.nf"

// Calculate hg_size from region param
def calc_hg_size(region) {
    def matcher = region =~ /[^:]+:(\d+)-(\d+)/
    if (!matcher) return ''
    def size = matcher[0][2].toLong() - matcher[0][1].toLong()
    if (size >= 1_000_000_000) return "${Math.round(size / 1_000_000_000)}g"
    if (size >= 1_000_000)     return "${Math.round(size / 1_000_000)}m"
    return "${Math.max(1, Math.round(size / 1000))}k"
}






// Define entry workflow and channels
workflow {
    def hg_size = calc_hg_size(params.region)
// This initial step is reading in the rows of the metadata file, splitting
    samples_ch = Channel
    .fromPath(params.metadata_table)
    .splitCsv(sep: '\t', header: ['sample_id', 'sex'], skip: 2)
    .map { row -> tuple(row.sample_id.toString().trim(), row.sex.toString().trim()) }


// Channel for finding the input aligned bam file
// Resolve globs explicitly to avoid empty/ambiguous file() behavior
    input_bam_ch = samples_ch
        .map { sample_id, sex ->
            def bam_pattern = "${params.bam_dir}/*/${sample_id}*${params.input_suffix}"
            def bai_pattern = "${params.bam_dir}/*/${sample_id}*${params.input_suffix}.bai"
            def bam_matches = files(bam_pattern)
            def bai_matches = files(bai_pattern)

            if (!bam_matches || bam_matches.size() == 0) {
                error "No BAM matched for ${sample_id} (${sex}) with pattern: ${bam_pattern}"
            }
            if (!bai_matches || bai_matches.size() == 0) {
                error "No BAI matched for ${sample_id} (${sex}) with pattern: ${bai_pattern}"
            }
            if (bam_matches.size() > 1) {
                error "Multiple BAM matches for ${sample_id} (${sex}): ${bam_matches}"
            }
            if (bai_matches.size() > 1) {
                error "Multiple BAI matches for ${sample_id} (${sex}): ${bai_matches}"
            }

            tuple(sample_id, sex, bam_matches[0], bai_matches[0])
        }


// Extract reads from the input aligned bam file
    extract_reads(input_bam_ch)

    // Branch based on sex
    extract_reads.out.fastq
        .map { sample_id, sex, fastq -> tuple(sample_id, sex, fastq, hg_size) }
        .filter { sample_id, sex, fastq, size -> sex == "XY" }
        .set { xy_fastq_ch }

    extract_reads.out.fastq
        .map { sample_id, sex, fastq -> tuple(sample_id, sex, fastq, hg_size) }
        .filter { sample_id, sex, fastq, size -> sex == "XX" }
        .set { xx_fastq_ch }
    

    run_hifiasm_XY(xy_fastq_ch)
    run_hifiasm_XX(xx_fastq_ch)


    // Combine the hap1 and hap2 fasta file into a channel to then be concatenated
    run_hifiasm_XX.out.hap1_fa
        .join(run_hifiasm_XX.out.hap2_fa, by: [0, 1])
        .set { xx_assemblies_ch}
    
    // concatenate hap1 and hap2 from run_hifiasm_XX
    cat_assemblies(xx_assemblies_ch)

    // add the fastq file to the tuple from extract reads for XX samples
    cat_assemblies.out.combined_diploid_fa
        .join(extract_reads.out.fastq, by: [0, 1])
        .set { xx_assembly_and_fastq_for_alignment_ch }
    
    // add the fastq file to the tuple from extract reads for XY samples
    run_hifiasm_XY.out.primary_fa
        .join(extract_reads.out.fastq, by: [0, 1])
        .set { xy_assembly_and_fastq_for_alignmet_ch }

    // combine XX and XY samples into one channel for alignment
    // new tuple will be sample_id, sex, assembly_fa, fastq
    xx_assembly_and_fastq_for_alignment_ch
        .mix(xy_assembly_and_fastq_for_alignmet_ch)
        .set { joined_assemblies_and_fastqs_for_alignment_ch }


    // outputs are (sample_id, sex, bam_path) as reads_to_assembly_bam
    // and (sample_id, sex, bam_index) as reads_to_assembly_bai
    align_to_assembly(joined_assemblies_and_fastqs_for_alignment_ch)

    // set up channels for exonerate
    // add XY primary hap label
    run_hifiasm_XY.out.primary_fa
        .map { sample_id, sex, primary_fa -> tuple(sample_id, sex, primary_fa, "primary") }
        .set { xy_assembly_for_exonerate_ch }
    
    run_hifiasm_XX.out.hap1_fa
        .map { sample_id, sex, hap1_fa -> tuple(sample_id, sex, hap1_fa, "hap1") }
        .set { xx_hap1_for_exonerate_ch }

    run_hifiasm_XX.out.hap2_fa
        .map { sample_id, sex, hap2_fa -> tuple(sample_id, sex, hap2_fa, "hap2") }
        .set { xx_hap2_for_exonerate_ch }

    // tuple structure is now sample_id, sex, assembly_fa, hap_name
    xy_assembly_for_exonerate_ch
        .mix(xx_hap1_for_exonerate_ch)
        .mix(xx_hap2_for_exonerate_ch)
        .set { all_assemblies_for_exonerate_ch }

    // Run exonerate on tuple of sample_id, sex, assembly_fa and hap_name
    // Output a tupe of sample_id, sex, hap_name, and exonerate_output_vulgar
    run_exonerate(all_assemblies_for_exonerate_ch)

    // run convert_vulgar_to_gff and output a tuple of sample_id, sex, hap_name and the vulgar_to_gff_output
    convert_vulgar_to_gff(run_exonerate.out.exonerate_output_vulgar)

    // run convert gff to bed and output a tuple of sample_id, sex, hap_name, and the gff_to_bed_output
    convert_gff_to_bed(convert_vulgar_to_gff.out.vulgar_to_gff_output)

    // Create new channel with sample_id, sex, hap_name, gff_to_bed_output, reads_to_assembly_bam, reads_to_assembly_bai
    convert_gff_to_bed.out.gff_to_bed_output
        .combine(align_to_assembly.out.reads_to_assembly_bam, by: [0,1])
        .combine(align_to_assembly.out.reads_to_assembly_bai, by: [0,1])
        .set { beds_and_bams_for_adding_stats_ch }
    
    // add alignment stats and output a tuple of sample_id, sex, hap_name, and stats_gff
    add_alignment_stats(beds_and_bams_for_adding_stats_ch)

    analyze_haplotype(add_alignment_stats.out.stats_gff)

    analyze_haplotype.out.haplotype_analysis
        .map { sample_id, sex, hap_name, tsv -> tsv }
        .collect()
        .set { all_tsvs_ch }
    
    concatenate_results(all_tsvs_ch)

    // parse gffs from exonerate to take coordinates of first 2 genes
    // output tuple containing sample_id, sex, hap_name, and a gene_coordinate_tsv
    parse_gff(convert_vulgar_to_gff.out.vulgar_to_gff_output)


 // rearrange the assemblies channel tuple so we can combine by sample_id, sex, and hap_name
    all_assemblies_for_exonerate_ch
        .map { sample_id, sex, assembly_fa, hap_name -> tuple(sample_id, sex, hap_name, assembly_fa) }
        .set { all_assemblies_by_hap_ch }


    // take the output tuple from parse_gff and map the entries from the coordinate tsv to it
    parse_gff.out.gene_coordinate_tsv
        .flatMap { sample_id, sex, hap_name, tsv -> 
            tsv.splitCsv(sep: '\t', header: true)
                .collect { row ->
                    tuple(sample_id, sex, hap_name, row.contig, row.start_position, row.end_position, row.type, row.gene_rank)
                }
        }
        .set { gene_rows_ch }
    

    // combine the gene rows with the assemblies by sample_id, sex, hap_name
    gene_rows_ch
        .combine(all_assemblies_by_hap_ch, by: [0, 1, 2])
        .set { gene_rows_with_assemblies_ch }



    // This process outputs a tuple that contains sample_id, sex, hap_name, contig, type, gene_rank, and the extracted fasta
    extract_annotation_sequences(gene_rows_with_assemblies_ch)


    // split the outputs from extracting fastas into gene 1 and 2
    extract_annotation_sequences.out.extracted_fasta
        .filter { sample_id, sex, hap_name, contig, type, gene_rank, extracted_fasta -> gene_rank == "gene1" }
        .map { sample_id, sex, hap_name, contig, type, gene_rank, extracted_fasta -> 
            tuple(sample_id, sex, hap_name, contig, type, gene_rank, file(params.OPN1LW_reference), extracted_fasta) }
        .set { gene1_extracted_annotation_fasta_with_ref_ch}
    
    extract_annotation_sequences.out.extracted_fasta
        .filter { sample_id, sex, hap_name, contig, type, gene_rank, extracted_fasta -> gene_rank == "gene2" }
        .map { sample_id, sex, hap_name, contig, type, gene_rank, extracted_fasta ->
            tuple(sample_id, sex, hap_name, contig, type, gene_rank, file(params.OPN1MW_reference), extracted_fasta) }
        .set { gene2_extracted_annotation_fasta_with_ref_ch}

    // add the path to the OPN1LW reference for gene1 dipcall runs and OPN1MW for gene2 dipcall runs
    // dipcall will expect a tuple of sample_id, sex, hap_name, contig, type, 
    gene1_extracted_annotation_fasta_with_ref_ch
        .mix(gene2_extracted_annotation_fasta_with_ref_ch)
        .set { all_genes_for_dipcall_ch }

    run_dipcall(all_genes_for_dipcall_ch)

    convert_vcf_coordinates(run_dipcall.out.dipcall_pair_vcf)

    run_vep(convert_vcf_coordinates.out.converted_dipcall_vcf)

    // run vep 115 on the converted-coordinate vcf and output a tuple with:
    // sample_id, sex, hap_name, contig, type, gene_rank, vep115_raw_vcf

    // Filter the output typle for gene 1 and 2
    // Gene 1 tuples have OPN1LW added and Gene 2 tuples have OPN1MW added
    run_vep.out.vep115_raw_vcf
        .filter { sample_id, sex, hap_name, contig, type, gene_rank, vep115_raw_vcf -> gene_rank == "gene1" }
        .map { sample_id, sex, hap_name, contig, type, gene_rank, vep115_raw_vcf ->
            tuple(sample_id, sex, hap_name, contig, type, gene_rank, "OPN1LW", vep115_raw_vcf) }
            .set { gene1_vep_outputs_L_ref_ch }
    
    run_vep.out.vep115_raw_vcf
        .filter { sample_id, sex, hap_name, contig, type, gene_rank, vep115_raw_vcf -> gene_rank == "gene2" }
        .map { sample_id, sex, hap_name, contig, type, gene_rank, vep115_raw_vcf ->
            tuple(sample_id, sex, hap_name, contig, type, gene_rank, "OPN1MW", vep115_raw_vcf) }
            .set { gene2_vep_outputs_M_ref_ch }   
    
    // Combine gene 1 and gene 2 tuples
    // New tuple structure is sample_id, sex, hap_name, contig, type, gene_rank, (gene_ref) OPN1LW|OPN1MW, vep115_raw_vcf
    gene1_vep_outputs_L_ref_ch
        .mix(gene2_vep_outputs_M_ref_ch)
        .set { vep_outputs_with_ref_genes_for_processing_ch }

    
    // run process vep vcf on tuple of: 
    // sample_id, sex, hap_name, contig, type, gene_rank, (gene_ref) OPN1LW|OPN1MW, vep115_raw_vcf
    // output new tuple as sample_id, sex, hap_name, contig, type, gene_rank, gene_ref, processed_vep_tsv
    process_vep_vcf(vep_outputs_with_ref_genes_for_processing_ch)

    // condense output tuple from process vep_vcf to
    // processed_vep_tsv
    process_vep_vcf.out.processed_vep_tsv
        .map { sample_id, sex, hap_name, contig, type, gene_rank, gene_ref, processed_vep_tsv -> processed_vep_tsv }
        .collect()
        .set { all_processed_vep_tsvs_ch }
    
    combine_SNP_analysis(all_processed_vep_tsvs_ch)

}
