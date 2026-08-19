# Introduction
This repository contains a workflow that resolves the gene copy-number, order, phasing, and variant calling for the opsin genes located at chromosome Xq28. The assembly and annotation steps of this workflow were used in (Anderson et al., 2026) Long-read sequencing with targeted assembly of the opsin locus accurately evaluates genes in expressed positions. https://www.medrxiv.org/content/10.64898/2026.03.17.26348636v1

## Usage
This workflow is written for nextflow version 26
```
nextflow run main.nf \
--bam_dir \  path to directory of bam files
--input_suffix \  input bam file suffix
--region_name \  name or identifier for samples in batch
--region \  coordinates in chr-start-end format 
--metadata_table \  list of samples
--output_dir \  publish location for outputs
--final_output_name \  final summary file name
--nested_bams \  flag that looks for bams in nested sub-directories ${bam_dir}/*/${sample_id}*${input_suffix} (default is FALSE; to run just add the flag to the nextflow command)
-resume
```

### Note about input bam directory
- If the bams are nested, use the --nested_bams flag
## Dependencies and environments
All dependencies should be available through conda, docker or github
- samtools - version 1.22 or newer
- hifiasm
- minimap2 - version 2.28 or newer
- exonerate
- vep - version 115
- dipcall - version 0.3 (follow instructions from github)


## dipcall edits for usage in workflow:
- For dipcall to work properly, you must open the dipcall-aux.js file and change line 160 to (min_var_len  = 10000)
## Input and output file formats
The starting file for this workflow must be an aligned bam file. The bam index is also needed. Note that the reference genome will effect the genomic coordinates that you use. 

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

2. A summary SNV file where each sample, haplotype and first two annotated genes are appended to. This file will be named using the final output name and will end in "combined_SNP_analysis.tsv".

|Column name|Contents|
|-----------|--------|
|sample|Sample identifier from the original metadata file|
|sex|Sex of the sample (can be XX or XY)|
|haplotype|hap1 or hap2 for XX and primary for XY samples|
|gene_rank|The order of first two genes annotated (gene1 or gene2)|
|gene_ref|The reference gene that the annotated gene has its variants called against (gene1 -> OPN1LW gene2 -> OPN1MW)|
|gene_annotation|The gene that was annotated in the contig (OPN1LW_exon5 or OPN1LW_exon5) This will tell you if you have an L or M annotation in the first or second position|
|Codons 65-309|These columns will have the reference nucleotides for the codon in the gene_ref. If there are no variants, all letters are capitalized (AGA). If there is a variant then the capitalized letter will be the variant (AGA -> agG with A->G being the variant)|
|AA|This is the translation of all the amino acids from the codon list|
|exon3_combo|Combination of codons 153, 171, 174, 178, and 180 in exon 3|
