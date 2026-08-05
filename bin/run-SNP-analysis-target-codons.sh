#!/bin/bash

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <sample_list> <vep_output_dir> <snp_analysis_output_dir> <summary_file_name>"
    exit 1
fi

SAMPLE_LIST="$1"
VEP_OUTPUT_DIR="$2"
SNP_ANALYSIS_OUTPUT_DIR="$3"
SUMMARY_FILE_NAME="$4"
scripts_dir="/n/zanderson/OPSIN-carrier-screen/manual-assembly-and-annotate/scripts"
mkdir -p "$SNP_ANALYSIS_OUTPUT_DIR"

process_vep_vcf() {
    local vep_vcf="$1"
    local gene_ref="$2"

    python3 $scripts_dir/opsin_vep_snp_analysis_target_codons.py "$vep_vcf" "$gene_ref" "$SNP_ANALYSIS_OUTPUT_DIR"
}

while IFS=$'\t' read -r SAMPLE_ID SEX EXTRA_COLS; do
    [[ "$SAMPLE_ID" == "SampleID" || "$SAMPLE_ID" == \#* || -z "$SAMPLE_ID" ]] && continue ## skip header and empty lines

    VEP_SAMPLE_DIR="${VEP_OUTPUT_DIR}/${SAMPLE_ID}_${SEX}_VEP"
    if [[ "$SEX" == "XX" ]]; then
        for HAP_NAME in hap1 hap2; do
            for GENE_RANK in gene1 gene2; do
                for TYPE in OPN1LW_exon5 OPN1MW_exon5; do
                    vcf_file="${VEP_SAMPLE_DIR}/${SAMPLE_ID}_${SEX}_${HAP_NAME}_${GENE_RANK}_${TYPE}"*".vep.114.vcf"
                    for f in $vcf_file; do
                        [ -e "$f" ] || continue
                        if [ "$GENE_RANK" == "gene1" ]; then
                        process_vep_vcf "$f" "OPN1LW"
                        elif [ "$GENE_RANK" == "gene2" ]; then
                            process_vep_vcf "$f" "OPN1MW"
                        fi
                    done
                done
            done
        done

    elif [[ "$SEX" == "XY" ]]; then
        for GENE_RANK in gene1 gene2; do
            for TYPE in OPN1LW_exon5 OPN1MW_exon5; do
                vcf_file="${VEP_SAMPLE_DIR}/${SAMPLE_ID}_${SEX}_primary_${GENE_RANK}_${TYPE}"*".vep.114.vcf"
                for f in $vcf_file; do
                    [ -e "$f" ] || continue
                    if [ "$GENE_RANK" == "gene1" ]; then
                        process_vep_vcf "$f" "OPN1LW"
                    elif [ "$GENE_RANK" == "gene2" ]; then
                        process_vep_vcf "$f" "OPN1MW"
                    fi
                done
            done
        done
    fi

done < "$SAMPLE_LIST"

# Combine all individual TSVs into one
python3 $scripts_dir/combine_SNP_analysis2.py "$SNP_ANALYSIS_OUTPUT_DIR" "${SNP_ANALYSIS_OUTPUT_DIR}/${SUMMARY_FILE_NAME}"