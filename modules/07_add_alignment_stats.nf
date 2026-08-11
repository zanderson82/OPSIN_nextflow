process add_alignment_stats {

    tag "${sample_id}-${sex}-${params.region_name}-${hap_name}"
    label "add_alignment_stats"

    input:
    tuple val(sample_id), val(sex), val(hap_name), path(gff_to_bed_output), path(reads_to_assembly_bam), path(reads_to_assembly_bai)

    output:
    tuple val(sample_id), val(sex), val(hap_name), path("${sample_id}-${sex}-${params.region_name}-${hap_name}_stats.gff"), emit: stats_gff

    script:
    """
    echo -e "contig\\tstart\\tend\\ttype\\tstrand\\treads\\tMQ0\\tratio" > "${sample_id}-${sex}-${params.region_name}-${hap_name}_stats.gff"

    if [ ! -s "${gff_to_bed_output}" ]; then
        exit 0
    fi

    tail -n +2 "${gff_to_bed_output}" | while IFS=\$'\\t' read -r contig start end type score strand; do
        [ -z "\$contig" ] && continue

        safe_type=\$(echo "\$type" | tr '/' '_')
        region_prefix="temp_\${contig}-\${start}-\${end}-\${safe_type}"
        samtools_start=\$((start + 1))

        samtools view -b -F 0x900 "${reads_to_assembly_bam}" "\${contig}:\${samtools_start}-\${end}" > "\${region_prefix}.bam" 2>/dev/null || true

        if [ -s "\${region_prefix}.bam" ]; then
            samtools index "\${region_prefix}.bam" 2>/dev/null || true
            samtools stats "\${region_prefix}.bam" > "\${region_prefix}.bam.stats" 2>/dev/null || true

            reads=\$(awk -F '\\t' '/^SN/ && /reads mapped:/{print \$3}' "\${region_prefix}.bam.stats" | head -n 1)
            mq0=\$(awk -F '\\t' '/^SN/ && /reads MQ0:/{print \$3}' "\${region_prefix}.bam.stats" | head -n 1)

            if [ -n "\$reads" ] && [ "\$reads" -ne 0 ] && [ -n "\$mq0" ]; then
                ratio=\$(awk -v m="\$mq0" -v r="\$reads" 'BEGIN{printf "%.10f", m/r}')
            else
                reads="\${reads:-0}"
                mq0="\${mq0:-0}"
                ratio="NA"
            fi
        else
            reads="0"
            mq0="0"
            ratio="NA"
        fi

        rm -f "\${region_prefix}.bam" "\${region_prefix}.bam.bai" "\${region_prefix}.bam.stats"
        echo -e "\${contig}\\t\${start}\\t\${end}\\t\${type}\\t\${strand}\\t\${reads}\\t\${mq0}\\t\${ratio}" >> "${sample_id}-${sex}-${params.region_name}-${hap_name}_stats.gff"
    done
    """
}