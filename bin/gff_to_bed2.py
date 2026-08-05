#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import re
import os

def write_empty_bed_outputs(sample_name, hap_name, output_dir):
    """Emit expected BED outputs when the GFF has no annotations."""
    os.makedirs(output_dir, exist_ok=True)
    main_output_file = os.path.join(output_dir, f"{sample_name}_{hap_name}.bed")
    with open(main_output_file, 'w') as f:
        f.write("contig\tstart\tend\ttype\tscore\tstrand\n")

    output_prefix = os.path.join(output_dir, f"{sample_name}_{hap_name}")
    open(f"{output_prefix}-LCR_and_exon5-annotation-file.bed", 'w').close()
    open(f"{output_prefix}-other-annotation-file.bed", 'w').close()
    print(f"No valid annotations found; wrote empty BED outputs for {sample_name} {hap_name}")
    return main_output_file

def gff_to_bed(gff_file, sample_name, sex, hap_name, output_dir):
    try:
        # Check if input file exists
        if not os.path.exists(gff_file):
            raise FileNotFoundError(f"GFF file not found: {gff_file}")
        
        # Read GFF file with proper column names
        df = pd.read_csv(gff_file, sep='\t', header=None, 
                        names=['seqname', 'source', 'similarity', 'start', 'end', 'score', 'strand', 'blank', 'attribute'])
        
        # Remove empty rows
        df = df.dropna(subset=['seqname'])
        
        if df.empty:
            return write_empty_bed_outputs(sample_name, hap_name, output_dir)
        
        # Convert data types
        df['start'] = pd.to_numeric(df['start'], errors='coerce')
        df['end'] = pd.to_numeric(df['end'], errors='coerce')
        df['score'] = pd.to_numeric(df['score'], errors='coerce')
        df['strand'] = df['strand'].astype(str)
        df['blank'] = df['blank'].astype(str)
        
        # Extract query name from attribute column
        df['query_name'] = df['attribute'].str.extract(r'query=([^;]+)')
        
        # Check if we successfully extracted query names
        if df['query_name'].isna().all():
            print(f"Warning: No query names found in attribute column for {gff_file}")
            df['query_name'] = 'unknown'

        # Handle strand-specific coordinate adjustments
        # Note: BED format uses 0-based coordinates, GFF uses 1-based
        # Convert GFF 1-based to BED 0-based by subtracting 1 from start
        df['bed_start'] = df['start'] - 1
        df['bed_end'] = df['end']
    
        mask_LCR_and_exon5 = df['query_name'].isin(['OPN1LW_exon5', 'OPN1MW_exon5', 'LCR'])
        
        output_bed_LCR_and_exon5 = df[mask_LCR_and_exon5][['seqname', 'bed_start', 'bed_end', 'query_name', 'score', 'strand']].copy()
        output_bed_other = df[~mask_LCR_and_exon5][['seqname', 'bed_start', 'bed_end', 'query_name', 'score', 'strand']].copy()
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Build output prefix using sample_name and hap_name
        output_prefix = os.path.join(output_dir, f"{sample_name}_{hap_name}")
        
        # Write main BED file (this is what the bash script expects)
        # Format: seqname, start, end, type, score, strand
        main_output_file = os.path.join(output_dir, f"{sample_name}_{hap_name}.bed")
        
        # Add header for the main BED file (bash script expects this)
        with open(main_output_file, 'w') as f:
            f.write("contig\tstart\tend\ttype\tscore\tstrand\n")
        output_bed_LCR_and_exon5.to_csv(main_output_file, sep='\t', header=False, index=False, mode='a')
        
        print(f"Successfully wrote main BED file: {main_output_file} ({len(output_bed_LCR_and_exon5)} records)")
        
        # Write LCR and exon5 specific file
        output_file_LCR_and_exon5 = f"{output_prefix}-LCR_and_exon5-annotation-file.bed"
        output_bed_LCR_and_exon5.to_csv(output_file_LCR_and_exon5, sep='\t', header=False, index=False)
        print(f"Successfully converted {len(output_bed_LCR_and_exon5)} records to {output_file_LCR_and_exon5}")
        
        # Write other annotations file
        output_file_other = f"{output_prefix}-other-annotation-file.bed"
        output_bed_other.to_csv(output_file_other, sep='\t', header=False, index=False)
        print(f"Successfully converted {len(output_bed_other)} records to {output_file_other}")
        
        return main_output_file
        
    except Exception as e:
        print(f"ERROR in gff_to_bed: {str(e)}")
        print(f"Failed to process file: {gff_file}")
        return write_empty_bed_outputs(sample_name, hap_name, output_dir)


def main():
    parser = argparse.ArgumentParser(description="Convert GFF to BED")
    parser.add_argument("--gff_file", type=str, required=True, help="Path to the GFF file")
    parser.add_argument("--sample_name", type=str, required=True, help="Name of the sample")
    parser.add_argument("--sex", type=str, required=True, help="Sex of the sample")
    parser.add_argument("--hap_name", type=str, required=True, help="Haplotype name")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to the output directory")
    args = parser.parse_args()

    gff_to_bed(args.gff_file, args.sample_name, args.sex, args.hap_name, args.output_dir)


if __name__ == "__main__":
    main()