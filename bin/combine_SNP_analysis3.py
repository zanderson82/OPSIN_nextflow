import sys, os, glob

# Sort order within each sample: hap1 before hap2 (or primary for XY), gene1 before gene2
HAP_ORDER  = {'hap1': 0, 'hap2': 1, 'primary': 0}
GENE_ORDER = {'gene1': 0, 'gene2': 1}


def row_sort_key(fields):
    hap  = fields[2] if len(fields) > 2 else ''
    gene = fields[3] if len(fields) > 3 else ''
    return (HAP_ORDER.get(hap, 99), GENE_ORDER.get(gene, 99))


def combine(input_dir, output_file):
    pattern = os.path.join(input_dir, '**', '*.vep_SNP_analysis.tsv')
    files = sorted(glob.glob(pattern, recursive=True))

    if not files:
        print(f"No files matching pattern in {input_dir}")
        sys.exit(1)

    header = None
    # Group data rows by (sample_id, sex)
    sample_groups = {}

    for f in files:
        with open(f) as fh:
            lines = fh.read().splitlines()
        if len(lines) < 2:
            continue
        if header is None:
            header = lines[0]
        row = lines[1]
        fields = row.split('\t')
        sample_id = fields[0] if len(fields) > 0 else 'unknown'
        sex       = fields[1] if len(fields) > 1 else 'unknown'
        key = (sample_id, sex)
        if key not in sample_groups:
            sample_groups[key] = []
        sample_groups[key].append(row)

    all_data_rows = []

    for (sample_id, sex) in sorted(sample_groups.keys()):
        rows = sorted(
            sample_groups[(sample_id, sex)],
            key=lambda r: row_sort_key(r.split('\t'))
        )

        # Per-sample combined file inside the existing sample/sex subdirectory
        per_sample_dir = os.path.join(input_dir, sample_id, sex)
        os.makedirs(per_sample_dir, exist_ok=True)
        per_sample_file = os.path.join(per_sample_dir, f'{sample_id}_{sex}_combined_SNP_analysis.tsv')
        with open(per_sample_file, 'w') as out:
            out.write(header + '\n')
            out.write('\n'.join(rows) + '\n')

        all_data_rows.extend(rows)

    # Global combined file
    with open(output_file, 'w') as out:
        out.write(header + '\n')
        out.write('\n'.join(all_data_rows) + '\n')

    print(f"Combined {len(files)} files -> {output_file}")
    print(f"Per-sample files written to {input_dir}/<sample>/<sex>/<sample>_<sex>_combined_SNP_analysis.tsv")


if len(sys.argv) != 3:
    print("Usage: python combine_SNP_analysis2.py <input_dir> <output_file>")
    print("Example: python combine_SNP_analysis2.py /path/to/snp_output combined_SNP_analysis.tsv")
    sys.exit(1)

combine(sys.argv[1], sys.argv[2])
