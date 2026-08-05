#!/usr/bin/env python3
import sys, re, os

TARGET_POSITIONS = {65, 100, 111, 116, 151, 153, 155, 171, 174, 178, 180, 203, 230, 233, 236, 277, 285, 309}
COMBINATION_POSITIONS = {153, 171, 174, 178, 180}

# Reference codons for each gene at each target protein position.
# Used when a position has no variant in the VCF (i.e. the sample matches the reference).
# Codons verified from symbol-matched VEP annotations; invariant positions derived from REF_AAS.
REF_NUCLEOTIDES = {
    'OPN1LW': {
        65:  'ACT',  # T
        100: 'CTA',  # L
        111: 'ATT',  # I
        116: 'TCT',  # S
        151: 'AGA',  # R
        153: 'ATG',  # M
        155: 'GTC',  # V
        171: 'GTG',  # V
        174: 'GCC',  # A
        178: 'ATC',  # I
        180: 'GCT',  # A
        203: 'TGC',  # C
        230: 'ATC',  # I
        233: 'GCT',  # A
        236: 'ATG',  # M
        277: 'TAC',  # Y
        285: 'ACC',  # T
        309: 'TAC',  # Y
    },
    'OPN1MW': {
        65:  'ATT',  # I
        100: 'CTG',  # L
        111: 'GTT',  # V
        116: 'TAT',  # Y
        151: 'AGA',  # R
        153: 'ATG',  # M
        155: 'GTC',  # V
        171: 'GTG',  # V
        174: 'GCC',  # A
        178: 'ATC',  # I
        180: 'GCT',  # A
        203: 'TGC',  # C
        230: 'ACC',  # T
        233: 'AGC',  # S
        236: 'GTG',  # V
        277: 'TTC',  # F
        285: 'GCC',  # A
        309: 'TTC',  # F
    },
}

REF_AAS = {
    'OPN1LW': {
        65: 'T',
        100: 'L',
        111: 'I',
        116: 'S',
        151: 'R',
        153: 'M',
        155: 'V',
        171: 'V',
        174: 'A',
        178: 'I',
        180: 'A',
        203: 'C',
        230: 'I',
        233: 'A',
        236: 'M',
        277: 'Y',
        285: 'T',
        309: 'Y',
    },
    'OPN1MW': {
        65: 'I',
        100: 'L',
        111: 'V',
        116: 'Y',
        151: 'R',
        153: 'M',
        155: 'V',
        171: 'V',
        174: 'A',
        178: 'I',
        180: 'A',
        203: 'C',
        230: 'T',
        233: 'S',
        236: 'V',
        277: 'F',
        285: 'A',
        309: 'F',
    },
}

SORTED_POSITIONS = sorted(TARGET_POSITIONS)
SORTED_COMBO = sorted(COMBINATION_POSITIONS)

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def merge_alt_codon(ref_codon_str, hit_list):
    """Merge multiple VEP alt codons that affect the same codon.

    VEP marks the changed base with uppercase in both ref and alt codon strings.
    Each variant in hit_list may change a different base within the codon.
    We apply all changes sequentially; the final string preserves uppercase at
    every position that was altered by any variant.

    ref_codon_str: pure-uppercase reference codon, e.g. 'GTG'
    hit_list:      list of hit dicts each carrying 'alt_codon' in VEP mixed-case
                   format, e.g. 'Atg' or 'gtT'
    Returns:       merged alt codon, e.g. 'AtT'
    """
    result = list(ref_codon_str.lower())
    for hit in hit_list:
        alt_vep = hit.get('alt_codon', '')
        if not alt_vep or len(alt_vep) != 3:
            continue
        for i, ch in enumerate(alt_vep):
            if ch.isupper():
                result[i] = ch
    return ''.join(result)


def parse_vcf_file(vcf_file, gene_ref):
    csq_fields = None
    hits = {}  # pos_int -> list of hit dicts (multiple variants can share a codon)

    if gene_ref not in REF_AAS:
        print(f"Error: unknown gene reference '{gene_ref}'. Choose from: {list(REF_AAS.keys())}")
        sys.exit(1)

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):
                csq_fields = re.search(r'Format: (.+?)"', line).group(1).split('|')
                continue
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            info = fields[7]
            if 'CSQ=' not in info:
                continue

            csq_str = info.split('CSQ=')[1].split(';')[0]
            vals = dict(zip(csq_fields, csq_str.split('|')))

            protein_position = vals.get('Protein_position', '')
            aa = vals.get('Amino_acids', '')
            consequence = vals.get('Consequence', '')
            gene = vals.get('SYMBOL', '')
            codons = vals.get('Codons', '')

            position_number = protein_position.split('/')[0]
            if not position_number.isdigit():
                continue

            pos_int = int(position_number)

            if pos_int in TARGET_POSITIONS or consequence == 'stop_gained':
                alt_aa = aa.split('/')[1] if '/' in aa else aa
                alt_codon = codons.split('/')[1] if '/' in codons else codons
                hit = {
                    'chrom': fields[0],
                    'pos': fields[1],
                    'ref': fields[3],
                    'alt': fields[4],
                    'gene': gene,
                    'consequence': consequence,
                    'amino_acids': aa,
                    'alt_aa': alt_aa,
                    'alt_codon': alt_codon,
                    'stop_gained': consequence == 'stop_gained',
                }
                hits.setdefault(pos_int, []).append(hit)

    ref_aas = REF_AAS[gene_ref]
    ref_nucs = REF_NUCLEOTIDES[gene_ref]
    result_aa = {}
    result_codon = {}

    for pos in SORTED_POSITIONS:
        ref_aa = ref_aas.get(pos, '?')
        ref_codon = ref_nucs.get(pos, '???')
        if pos in hits:
            hit_list = hits[pos]
            if any(h['stop_gained'] for h in hit_list):
                result_aa[pos] = '*'
                stop_hit = next(h for h in hit_list if h['stop_gained'])
                result_codon[pos] = stop_hit['alt_codon'] if stop_hit['alt_codon'] else ref_codon
            elif len(hit_list) == 1:
                hit = hit_list[0]
                result_aa[pos] = hit['alt_aa']
                result_codon[pos] = hit['alt_codon'] if hit['alt_codon'] else ref_codon
            else:
                # Multiple variants affect the same codon — merge all changes
                merged = merge_alt_codon(ref_codon, hit_list)
                result_codon[pos] = merged
                result_aa[pos] = CODON_TABLE.get(merged.upper(), '?')
        else:
            result_aa[pos] = ref_aa
            result_codon[pos] = ref_codon

    return hits, result_aa, result_codon


def write_output(vcf_file, gene_ref, output_file):
    # Derive sample name from filename:
    # e.g. GM19038_XX_hap1_gene1_OPN1LW_exon5_h1tg000001l.vep.114.vcf
    #   -> GM19038_XX_hap1_gene1_OPN1LW_exon5
    basename = os.path.basename(vcf_file)
    name = basename.replace('.vep.114.vcf', '').replace('.vep.vcf', '')
    sample_name = name

    # Parse metadata from filename: SAMPLEID_SEX_HAP_GENERANK_GENEREF_EXON
    name_parts = sample_name.split('-')
    sample_id       = name_parts[0] if len(name_parts) > 0 else 'unknown'
    sex             = name_parts[1] if len(name_parts) > 1 else 'unknown'
    haplotype       = name_parts[2] if len(name_parts) > 2 else 'unknown'
    gene_rank       = name_parts[3] if len(name_parts) > 3 else 'unknown'
    gene_annotation = name_parts[4] if len(name_parts) > 4 else 'unknown'
    hits, result_aa, result_codon = parse_vcf_file(vcf_file, gene_ref)

    # exon3_combo is based on amino acids at the 5 combination positions
    combo = '-'.join(result_aa[p] for p in SORTED_COMBO)
    # AA column: full 18-position amino acid string
    aa_sequence = ''.join(result_aa[p] for p in SORTED_POSITIONS)

    pos_cols = [str(p) for p in SORTED_POSITIONS]
    header = '\t'.join(
        ['sample', 'sex', 'haplotype', 'gene_rank', 'gene_ref', 'gene_annotation']
        + pos_cols
        + ['AA', 'exon3_combo']
    )

    codon_vals = [result_codon[p] for p in SORTED_POSITIONS]
    row = '\t'.join(
        [sample_id, sex, haplotype, gene_rank, gene_ref, gene_annotation]
        + codon_vals
        + [aa_sequence, combo]
    )

    out_path = output_file
    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(out_path, 'w') as f:
        f.write(header + '\n')
        f.write(row + '\n')

    print(f"Written: {out_path}")


if len(sys.argv) != 4:
    print("Usage: python vcf-SNP-analysis3.py <vcf_file> <gene_ref> <output_file>")
    print("Example: python vcf-SNP-analysis3.py sample.vep.114.vcf OPN1LW sample.vep_SNP_analysis.tsv")
    sys.exit(1)

write_output(sys.argv[1], sys.argv[2], sys.argv[3])
