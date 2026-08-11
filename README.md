# Introduction
This repository contains a workflow that resolves the gene copy-number, order, phasing, and variant calling for the opsin genes located at chromosome Xq28. The assembly and annotation steps of this workflow were used in (Anderson et al., 2026) Long-read sequencing with targeted assembly of the opsin locus accurately evaluates genes in expressed positions. https://www.medrxiv.org/content/10.64898/2026.03.17.26348636v1

## Usage
This workflow is written for nextflow version 26
```
nextflow run main.nf \
--bam_dir "$path_to_directory_of_bam_files \
--input_suffix "$input_bam_file_suffix" \
--region_name "$name_or_identifier_for_samples_in_batch"\
--region "$chr-start-end" \
--metadata_table "$list_of_samples" \
--final_output_name "$final_summary_file_name" \
-resume
```

## Dependencies and environments
All dependencies should be available through conda and or docker
- samtools - version 1.22 or newer
- hifiasm
- minimap2 - version 2.28 or newer
- exonerate
- vep - version 115
- dipcall - version 0.3

## Input and output file formats
The starting file for this workflow must be an aligned bam file. Note that the reference genome will effect the genomic coordinates that you use. 

There are two summary output files:
1. A summary annotation file that each sample and haplotype is appended to. This file has the following columns: 

|Column name|Contents|
|-----------|--------|
|sample_id|Sample identifier from the original metadata file|
|sex|Sex of the sample (can be XX or XY)|
|haplotype|hap1 or hap2 for XX and primary for XY samples|
|structure|Order of genes annotated on a haplotype (e.g., L-M)|
|lw_count|Number of OPN1LW genes annotated on the haplotype|
|mw_count|Number of OPN1MW genes annotated on the haplotype|
|total_genes|Total number of genes annotated on the haplotype|
|lcr_count|Number of locus control region(s) (LCR) annotated on the haplotype|
|total_contigs|Number of contigs assigned to said haplotype with annotations|
|contigs_with_lcr|Number of contigs with LCR annotations|
|contigs_without_lcr|Number of contigs without an LCR annotation|
|orphan_genes|Genes annotated on contigs that don't have an LCR annotation|
|arrays_found|Number of LCR + L or M genes found|
|is_reverse|Indicator if the array was assembled in reverse|
|orientation_ambiguous|Indicates that there were an equal number of annotations on + and - strands|
|primary_contig|Name of contig that is marked as primary (LCR annotation + most gene annotations if there are multiple contigs with arrays)|
|primary_lcr_position|Coordinate of primary LCR annotation on its contig (most helpful if there are multiple LCR annotations)|
|primary_lcr_ratio|Ratio of mapq0 that map to the LCR annotation site to the total number of reads|
|primary_lcr_reads|Number of reads that map to the primary LCR annotation site|
|primary_lcr_mapq0|Number of reads that map to the primary LCR annotation site with a mapq score of 0|

