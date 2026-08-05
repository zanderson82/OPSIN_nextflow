#!/usr/bin/env python3
import re
import sys
import os
from collections import OrderedDict, defaultdict

# Optional globals for writing first two genes to TSV
GENES_TSV_PATH = None
SAMPLE_ID = ""
HAP_NAME = ""

def analyze_opsin_exon5(gff_file, sample_id, hap_name, sex, genes_tsv_path):
    # Initialize counters for exon5
    lw_exon5_count = 0
    mw_exon5_count = 0
    
    # Store annotations grouped by contig
    contigs = defaultdict(list)
    
    # Debug all contents of GFF file
    print(f"Analyzing GFF file: {gff_file}")
    
    # Parse the GFF file to find the annotations
    try:
        with open(gff_file, 'r') as f:
            content = f.read()
            print(f"File loaded, total characters: {len(content)}")
            
            # If file is empty
            if not content.strip():
                print("WARNING: GFF file is empty!")
                return 0, 0, "Unknown (Empty GFF)", 0
            
            lines = content.splitlines()
            print(f"Total lines in file: {len(lines)}")
            
            # reading the gff file line by line
            for line_count, line in enumerate(lines, 1):
                if line.startswith('#') or line.strip() == '':
                    continue
                
                # Print the raw line for debugging
                print(f"Line {line_count}: {line.strip()}")
                
                # split the line into fields
                fields = line.strip().split('\t')
                # check if the line has fewer than 9 fields
                if len(fields) < 9:
                    print(f"  Warning: Line {line_count} has fewer than 9 fields ({len(fields)} fields)")
                    continue
                
                # Extract information from GFF format
                contig_name = fields[0]  # Contig/target name
                start_pos = int(fields[3]) if fields[3].isdigit() else 0
                end_pos = int(fields[4]) if fields[4].isdigit() else 0
                strand = fields[6]
                attributes = fields[8]
                
                # Extract query name from attributes
                query_name = ""
                if "query=" in attributes:
                    query_match = re.search(r'query=([^;]+)', attributes)
                    if query_match:
                        query_name = query_match.group(1)
                
                # Debug info
                print(f"  Contig: {contig_name}, Query: {query_name}, Position: {start_pos}-{end_pos}, Strand: {strand}")
                
                # Only process OPSIN-related annotations
                if query_name in ['LCR', 'OPN1LW_exon5', 'OPN1MW_exon5']:
                    # Create annotation object
                    annotation = {
                        'type': query_name,
                        'position': start_pos,
                        'end_position': end_pos,
                        'strand': strand,
                        'contig': contig_name
                    }
                    
                    # Add to contig-specific list
                    contigs[contig_name].append(annotation)
                    
                    # Count exon5 annotations for totals
                    if query_name == 'OPN1LW_exon5':
                        lw_exon5_count += 1
                        print(f"    Found LW exon5 on contig {contig_name} at position {start_pos}")
                    elif query_name == 'OPN1MW_exon5':
                        mw_exon5_count += 1
                        print(f"    Found MW exon5 on contig {contig_name} at position {start_pos}")
                    elif query_name == 'LCR':
                        print(f"    Found LCR on contig {contig_name} at position {start_pos}")
                
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        return 0, 0, "Error", 0
    
    # Calculate total LCR count BEFORE printing results
    total_lcr_count = 0
    for contig_name, annotations in contigs.items():    # Iterate through each contig and its annotations
        lcr_annotations = [a for a in annotations if a['type'] == 'LCR']    # Extract LCR annotations from the contig
        total_lcr_count += len(lcr_annotations)    # Add the number of LCR annotations to the total count
    
    # Print results
    print(f"\n=== Sample {sample_id} ({hap_name}) OPSIN Gene Analysis ===")
    
    print(f"  OPN1LW: {lw_exon5_count} copies (based on exon 5 count)")
    print(f"  OPN1MW: {mw_exon5_count} copies (based on exon 5 count)")
    print(f"  Total OPSIN genes: {lw_exon5_count + mw_exon5_count} copies")
    print(f"  Total LCR count: {total_lcr_count} LCRs")
    
    # Array structure analysis by contig
    print("\n=== Array Structure Analysis by Contig ===")
    
    # Check if we have any annotations
    if not contigs:
        print("WARNING: No OPSIN annotations found in the GFF file!")
        return 0, 0, "Unknown (No annotations)", 0
    
    print(f"Found OPSIN annotations on {len(contigs)} contig(s):")
    for contig_name, annotations in contigs.items():
        lcr_count = len([a for a in annotations if a['type'] == 'LCR'])
        exon5_count = len([a for a in annotations if a['type'] in ['OPN1LW_exon5', 'OPN1MW_exon5']])
        print(f"  Contig {contig_name}: {len(annotations)} annotations ({lcr_count} LCRs, {exon5_count} exon5s)")
    
    # Separate contigs with and without LCRs
    lcr_contigs = []
    non_lcr_contigs = []
    
    for contig_name, annotations in contigs.items():
        lcr_annotations = [a for a in annotations if a['type'] == 'LCR']
        exon5_annotations = [a for a in annotations if a['type'] in ['OPN1LW_exon5', 'OPN1MW_exon5']]
        # Split contigs into LCR and non-LCR contigs and store in separate lists with their exon5 counts and annotations
        if lcr_annotations:
            lcr_contigs.append({
                'name': contig_name,
                'annotations': annotations,
                'lcr_count': len(lcr_annotations),
                'exon5_count': len(exon5_annotations),
                'lcr_annotations': lcr_annotations,
                'exon5_annotations': exon5_annotations
            })
        elif exon5_annotations:  # Only add if it has exon5 annotations
            non_lcr_contigs.append({
                'name': contig_name,
                'annotations': annotations,
                'exon5_count': len(exon5_annotations),
                'exon5_annotations': exon5_annotations
            })
    
    print(f"\nContigs with LCRs: {len(lcr_contigs)}")
    print(f"Contigs with exon5s but no LCRs: {len(non_lcr_contigs)}")
    
    if not lcr_contigs:
        print("WARNING: No contigs with LCR annotations found!")
        if non_lcr_contigs:
            print("Falling back to non-LCR contigs for gene assignment...")
            all_exon5 = []
            for contig in non_lcr_contigs:
                all_exon5.extend(contig['exon5_annotations'])
            all_exon5.sort(key=lambda x: x['position'])

            for i, gene in enumerate(all_exon5):
                print(f"  {i+1}. {gene['type']} at position {gene['position']} on contig {gene['contig']}")

            if GENES_TSV_PATH and all_exon5:
                try:
                    header_needed = not os.path.exists(GENES_TSV_PATH) or os.path.getsize(GENES_TSV_PATH) == 0
                    with open(GENES_TSV_PATH, 'a') as out:
                        if header_needed:
                            out.write("sample_id\tsex\thap_name\tcontig\tstart_position\tend_position\ttype\tgene_rank\n")
                        for rank, gene in enumerate(all_exon5[:2], 1):
                            start_pos = max(1, gene['position'] - 12010)
                            end_pos = gene['end_position'] + 2393
                            out.write(f"{SAMPLE_ID}\t{SEX}\t{HAP_NAME}\t{gene['contig']}\t{start_pos}\t{end_pos}\t{gene['type']}\tgene{rank}\n")
                    print(f"Wrote {min(2, len(all_exon5))} gene annotation(s) from non-LCR contigs to: {GENES_TSV_PATH}")
                except Exception as e:
                    print(f"Failed to write gene annotations TSV: {e}")

            gene_parts = []
            for gene in all_exon5[:2]:
                if gene['type'] == 'OPN1LW_exon5':
                    gene_parts.append("L")
                elif gene['type'] == 'OPN1MW_exon5':
                    gene_parts.append("M")
            structure = "-".join(gene_parts) + " (no LCR)"
            print(f"\nFinal structure (no LCR): {structure}")
            return lw_exon5_count, mw_exon5_count, structure, total_lcr_count

        return lw_exon5_count, mw_exon5_count, "Unknown (No LCR found)", total_lcr_count
    
    # Sort LCR contigs by number of exon5s (descending) to prioritize the one with most genes
    lcr_contigs.sort(key=lambda x: x['exon5_count'], reverse=True)
    primary_contig = lcr_contigs[0]
    
    print(f"\nPrimary contig (most genes): {primary_contig['name']} with {primary_contig['exon5_count']} exon5s and {primary_contig['lcr_count']} LCRs")
    
    # Analyze primary contig structure
    print(f"\n--- Analyzing Primary Contig {primary_contig['name']} ---")
    primary_structure = analyze_contig_structure(primary_contig['name'], primary_contig['lcr_annotations'], primary_contig['exon5_annotations'], primary_contig['annotations'])
    
    # Collect additional genes from non-LCR contigs
    additional_genes = []
    if non_lcr_contigs:
        print(f"\n--- Processing Additional Contigs without LCRs ---")
        for contig in non_lcr_contigs:
            print(f"Contig {contig['name']}: {contig['exon5_count']} exon5 genes")
            
            # Sort genes by position for consistent ordering
            sorted_genes = sorted(contig['exon5_annotations'], key=lambda x: x['position'])
            additional_genes.extend(sorted_genes)
            
            for i, gene in enumerate(sorted_genes):
                print(f"  {i+1}. {gene['type']} at position {gene['position']}")
    
    # Build final combined structure
    if primary_structure and additional_genes:
        # Extract the gene portion from primary structure (remove orientation info)
        base_structure = primary_structure
        orientation_info = ""
        if "(reverse compliment)" in primary_structure:
            base_structure = primary_structure.replace(" (reverse compliment)", "")
            orientation_info = " (reverse compliment)"
        
        # Add additional genes
        additional_structure = ""
        for gene in additional_genes:
            if gene['type'] == 'OPN1LW_exon5':
                additional_structure += "-L"
            elif gene['type'] == 'OPN1MW_exon5':
                additional_structure += "-M"
        
        final_structure = base_structure + additional_structure + orientation_info
        
        print(f"\nCombined Array Structure:")
        print(f"  Primary contig: {primary_structure}")
        if additional_genes:
            print(f"  Additional genes: {additional_structure}")
        print(f"  Final structure: {final_structure}")

        # If only one gene on the primary contig but we have additional exon5 genes
        # append the first additional gene as gene2 to the TSV so we don't miss it
        if GENES_TSV_PATH and primary_contig['exon5_count'] == 1 and len(additional_genes) >= 1:
            try:
                add_gene = additional_genes[0]
                if add_gene['position'] - 12010 < 0:
                    start_position = 1
                else:
                    start_position = add_gene['position'] - 12010
                gene2 = {
                    'type': add_gene['type'],
                    'contig': add_gene['contig'],
                    'start_position': start_position,
                    'end_position': add_gene['end_position'] + 2393
                }
                with open(GENES_TSV_PATH, 'a') as out:
                    out.write(f"{SAMPLE_ID}\t{SEX}\t{HAP_NAME}\t{gene2['contig']}\t{gene2['start_position']}\t{gene2['end_position']}\t{gene2['type']}\tgene2\n")
                print(f"Appended second gene annotation from additional contigs to: {GENES_TSV_PATH}")
            except Exception as e:
                print(f"Failed to append second gene annotation from additional contigs: {e}")
        
        return lw_exon5_count, mw_exon5_count, final_structure, total_lcr_count
    elif primary_structure:
        print(f"\nFinal structure: {primary_structure}")
        return lw_exon5_count, mw_exon5_count, primary_structure, total_lcr_count
    else:
        print("No valid array structure could be determined")
        return lw_exon5_count, mw_exon5_count, "Unknown (No valid arrays found)", total_lcr_count

def analyze_contig_structure(contig_name, lcr_annotations, contig_exon5, all_annotations):
    """Analyze the array structure for a single contig"""
    
    if not lcr_annotations or not contig_exon5:
        return None
    
    # Determine the dominant strand orientation for this contig
    minus_strand_count = sum(1 for a in all_annotations if a['strand'] == '-')
    plus_strand_count = sum(1 for a in all_annotations if a['strand'] == '+')
    is_reverse_orientation = minus_strand_count > plus_strand_count
    
    print(f"Contig {contig_name} strand orientation: {'+' if not is_reverse_orientation else '-'} strand dominant " +
          f"({plus_strand_count} + strands, {minus_strand_count} - strands)")
    
    # Sort annotations on this contig by position
    if is_reverse_orientation:
        # For reverse orientation, sort in descending order
        all_annotations.sort(key=lambda x: x['position'], reverse=True)
        contig_exon5.sort(key=lambda x: x['position'], reverse=True)
        lcr_annotations.sort(key=lambda x: x['position'], reverse=True)
        print(f"Contig {contig_name} is in REVERSE orientation - sorting positions in descending order")
    else:
        # For forward orientation, sort in ascending order
        all_annotations.sort(key=lambda x: x['position'])
        contig_exon5.sort(key=lambda x: x['position'])
        lcr_annotations.sort(key=lambda x: x['position'])
        print(f"Contig {contig_name} is in FORWARD orientation - sorting positions in ascending order")
    
    # Debug: Print sorted annotations for this contig
    relevant_annotations = [a for a in all_annotations if a['type'] in ['LCR', 'OPN1LW_exon5', 'OPN1MW_exon5']]
    print(f"Sorted relevant annotations on contig {contig_name}:")
    for i, a in enumerate(relevant_annotations):
        print(f"  {i+1}. {a['type']} at position {a['position']} on strand {a['strand']}")
    
    # Handle multiple LCRs on this contig by grouping genes with their closest LCR
    if len(lcr_annotations) > 1:
        print(f"Multiple LCR annotations detected on contig {contig_name} ({len(lcr_annotations)}). Analyzing separate arrays:")
        
        # Group each exon5 with its closest LCR on this contig
        arrays = []
        for lcr in lcr_annotations:
            lcr_position = lcr['position']
            # Find exon5s closest to this LCR
            closest_exon5s = []
            for exon5 in contig_exon5:
                # Find which LCR this exon5 is closest to
                distances = [(abs(exon5['position'] - l['position']), l) for l in lcr_annotations]
                closest_lcr = min(distances, key=lambda x: x[0])[1]
                if closest_lcr['position'] == lcr_position:
                    closest_exon5s.append(exon5)
            
            if closest_exon5s:
                arrays.append({
                    'lcr': lcr,
                    'exon5s': closest_exon5s
                })
        
        # Print array information
        for i, array in enumerate(arrays):
            print(f"  Array {i+1} on contig {contig_name}: LCR at {array['lcr']['position']} with {len(array['exon5s'])} genes")
            for j, exon5 in enumerate(array['exon5s']):
                print(f"    Gene {j+1}: {exon5['type']} at position {exon5['position']}")
        
        # Use the array with the most genes
        if arrays:
            primary_array = max(arrays, key=lambda x: len(x['exon5s']))
            primary_exon5s = primary_array['exon5s']
            
            print(f"Using primary array on contig {contig_name} (most genes) for structure determination:")
            print(f"  LCR at position {primary_array['lcr']['position']}")
            print(f"  {len(primary_exon5s)} associated genes")
            
            return build_array_structure(primary_exon5s, is_reverse_orientation)
    else:
        # Single LCR on this contig - use original logic
        first_lcr = lcr_annotations[0]
        lcr_position = first_lcr['position']
        
        print(f"Single LCR found on contig {contig_name} at position {lcr_position}")
        
        # Find the closest exon5 to the LCR based on absolute distance
        closest_exon5 = min(contig_exon5, key=lambda x: abs(x['position'] - lcr_position))
        
        print(f"Closest exon5 to LCR on contig {contig_name}: {closest_exon5['type']} at position {closest_exon5['position']}")
        
        # Reorder contig_exon5 to put the closest one first
        reordered_exon5 = [closest_exon5] + [e for e in contig_exon5 if e != closest_exon5]
        
        # Print the reordered annotations
        print(f"Reordered exon5 annotations on contig {contig_name}:")
        for i, ann in enumerate(reordered_exon5):
            print(f"  {i+1}. {ann['type']} at position {ann['position']} on strand {ann['strand']}")
        
        return build_array_structure(reordered_exon5, is_reverse_orientation)
    
    return None

def build_array_structure(exon5_list, is_reverse_orientation):
    """Build the array structure string from a list of exon5 annotations"""
    
    if not exon5_list:
        return "Unknown (no exon5 annotations found)"
    
    # Build array structure
    full_structure = "[LCR]-"
    for a in exon5_list:
        if a['type'] == 'OPN1LW_exon5':
            full_structure += "L-"
        elif a['type'] == 'OPN1MW_exon5':
            full_structure += "M-"
    
    # Remove trailing dash
    full_structure = full_structure.rstrip('-')
    
    # Add orientation indicator
    orientation_indicator = " (reverse compliment)" if is_reverse_orientation else ""
    
    # Build gene structure with consistent dashes
    gene_parts = []
    for a in exon5_list:
        if a['type'] == 'OPN1LW_exon5':
            gene_parts.append("L")
        elif a['type'] == 'OPN1MW_exon5':
            gene_parts.append("M")
    
    # Join with dashes
    gene_structure = "-".join(gene_parts)
    
    print(f"  Full array structure: {full_structure}")
    first_gene = None
    second_gene = None
    if len(exon5_list) >= 1:
        first_gene = {
            'type': exon5_list[0]['type'],
            'contig': exon5_list[0]['contig'],
            'start_position': exon5_list[0]['position'] - 12010,
            'end_position': exon5_list[0]['end_position'] + 2393
        }
    if len(exon5_list) >= 2:
        second_gene = {
            'type': exon5_list[1]['type'],
            'contig': exon5_list[1]['contig'],
            'start_position': exon5_list[1]['position'] - 12010,
            'end_position': exon5_list[1]['end_position'] + 2393
        }
    
    # Optionally write the first two gene annotations to a TSV
    if GENES_TSV_PATH and (first_gene is not None):
        try:
            header_needed = not os.path.exists(GENES_TSV_PATH) or os.path.getsize(GENES_TSV_PATH) == 0
            with open(GENES_TSV_PATH, 'a') as out:
                if header_needed:
                    out.write("sample_id\tsex\thap_name\tcontig\tstart_position\tend_position\ttype\tgene_rank\n")
                out.write(f"{SAMPLE_ID}\t{SEX}\t{HAP_NAME}\t{first_gene['contig']}\t{first_gene['start_position']}\t{first_gene['end_position']}\t{first_gene['type']}\tgene1\n")
                if second_gene is not None:
                    out.write(f"{SAMPLE_ID}\t{SEX}\t{HAP_NAME}\t{second_gene['contig']}\t{second_gene['start_position']}\t{second_gene['end_position']}\t{second_gene['type']}\tgene2\n")
            print(f"Wrote first two gene annotations to: {GENES_TSV_PATH}")
        except Exception as e:
            print(f"Failed to write first two gene annotations TSV: {e}")
    
    return f"{gene_structure}{orientation_indicator}"
    
if __name__ == "__main__":
    gff_file = sys.argv[1]
    sample_id = sys.argv[2]
    hap_name = sys.argv[3]
    sex = sys.argv[4]
    GENES_TSV_PATH = sys.argv[5]
    SAMPLE_ID = sample_id
    HAP_NAME = hap_name
    SEX = sex
    analyze_opsin_exon5(gff_file, sample_id, hap_name, sex, GENES_TSV_PATH)
