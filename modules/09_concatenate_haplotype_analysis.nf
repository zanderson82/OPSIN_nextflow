process concatenate_results {
    publishDir "${params.output_dir}", mode: 'copy'
    label "concatenate_results"

    input:
    path(tsvs)

    output:
    path("${params.final_output_name}_summary.tsv"), emit: summary_tsv

    script:
    """
    python_opsin_processing_V6.py --header > "${params.final_output_name}_summary.tsv"
    cat ${tsvs} >> "${params.final_output_name}_summary.tsv"
    python - <<'PY'
    import csv
    summary_path = "${params.final_output_name}_summary.tsv"
    with open(summary_path, newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        rows = list(reader)

        with open(summary_path, 'w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(rows[0])
            writer.writerows(rows[1:])
    PY
    """
}