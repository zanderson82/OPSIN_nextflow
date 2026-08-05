#!/usr/bin/env python3
import sys


def convert_opsin_vcf_coordinates(
    sample_id,
    hap_name,
    sex,
    gene_rank,
    type,
    contig,
    vcf_file,
    output_vcf
):
    """
    Read a dipcall *.pair.vcf, shift coordinates to genomic space for OPSIN exon5,
    and write a new VCF with \"_converted.vcf\" suffix in the same dipcall directory.
    """


    with open(vcf_file, "r") as fin, open(output_vcf, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                # Preserve all original VCF headers so zero-variant files stay valid.
                fout.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                # Keep malformed lines out of output.
                continue

            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                # Skip lines where POS is not numeric.
                continue

            # Coordinate conversion for OPSIN exon5 on chrX
            if gene_rank == "gene1" and chrom == "chrX":
                pos += 154144242
            elif gene_rank == "gene2" and chrom == "chrX":
                pos += 154182595

            fields[1] = str(pos)
            fout.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    # Expect 9 arguments after the script name
    if len(sys.argv) != 9:
        sys.stderr.write(
            "Usage: convert-opsin-vcf-coordinates.py "
            "SAMPLE_ID HAP_NAME SEX GENE_RANK "
            "TYPE CONTIG INPUT_VCF OUTPUT_VCF\n"
        )
        sys.exit(1)

    (
        sample_id,
        hap_name,
        sex,
        gene_rank,
        type_,
        contig,
        vcf_file,
        output_vcf
    ) = sys.argv[1:]

    convert_opsin_vcf_coordinates(
        sample_id,
        hap_name,
        sex,
        gene_rank,
        type_,
        contig,
        vcf_file,
        output_vcf
    )