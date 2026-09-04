# Introduction
This repository contains a workflow that resolves the gene copy-number, order, phasing, and variant calling for the opsin genes located at chromosome Xq28. The assembly and annotation steps of this workflow were used in (Anderson et al., 2026) Long-read sequencing with targeted assembly of the opsin locus accurately evaluates genes in expressed positions. https://www.medrxiv.org/content/10.64898/2026.03.17.26348636v1

# Introduction
This repository contains a workflow that resolves the gene copy-number, order, phasing, and variant calling for the opsin genes located at chromosome Xq28. The assembly and annotation steps of this workflow were used in (Anderson et al., 2026) Long-read sequencing with targeted assembly of the opsin locus accurately evaluates genes in expressed positions. https://www.medrxiv.org/content/10.64898/2026.03.17.26348636v1


## Installation and setup
1. Clone this repository to your machine
2. Install hifiasm (requiring g++ and zlib)
    - git clone https://github.com/chhylp123/hifiasm into the ``./resources`` directory.
    - cd hifiasm && make
3. Install dipcall ()
    - wget https://github.com/lh3/dipcall/releases/download/v0.3/dipcall-0.3_x64-linux.tar.bz2
    - tar -jxf dipcall-0.3_x64-linux.tar.bz2
3. Create conda environments for each yaml file in ``./resources/environments``
    - ``conda env create -f environment.yml``
    - Each environment will have a name that is designated in the .yml file.
4. Update environment paths in ``./nextflow.config``
    - See **Using conda environments below**
5. Create and activate an environment with nextflow version 26
    - e.g., ``conda create -n nextflow_26 -c conda-forge -c bioconda nextflow=26.04.4``
6. Update ``./resources/sample_file.tsv`` with your sample IDs and the Sex
    - See **Input and output file formats**


### Using conda environments
This workflow is designed to use conda environments for the different steps or modules (found in ``./modules/``).
The conda environment that is to be used for each step is denoted in ``./nextflow.config``. Once you create your own environments with the YAML files, you will change the paths to each respective environment.

Example:
This is what the ``./nextflow.config`` currently looks like, starting at line 22.

```
process {
    withLabel: 'extract_reads' {
        conda = '/usr/share/millerlab/samtools-1.22'
        cpus = 10
    }

    withLabel: 'run_hifiasm_XY' {
        conda = '/usr/share/millerlab/hifiasm-0.25.0'
        cpus = 10
    }

    withLabel: 'run_hifiasm_XX' {
        conda = '/usr/share/millerlab/hifiasm-0.25.0'
        cpus = 10
    }

    withLabel: 'align_to_assembly' {
        conda = '/usr/share/millerlab/minimap-2.28'
        cpus = 10
    }

    withLabel: 'run_exonerate' {
        conda = '/home/zanderson/.conda/envs/exonerate-env'
    }

    withLabel: 'convert_gff_to_bed' {
        conda = '/home/zanderson/.conda/envs/exonerate-env'
    }

    withLabel: 'analyze_haplotype' {
        conda = '/home/zanderson/.conda/envs/exonerate-env'
    }
    withLabel: 'concatenate_results' {
        conda = '/home/zanderson/.conda/envs/exonerate-env'
    }
    withLabel: 'run_vep' {
        conda = '/usr/share/millerlab/vep-115.2'
        cpus = 10
    }
}
```


## Input and output file formats

### Aligned bam file
The starting file for this workflow must be an **aligned bam file** along with its bam index. Note that the reference genome will effect the genomic coordinates that you use. 

### Metadata file
The other file you will need is a tab-separated metadata file that has the sample identifier (Sample ID) and the Sex (XX or XY)

### Output files

A summary annotation file that each sample and haplotype is appended to will be generated and output to the ``output_dir`` location. This file has the following columns: 

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


## Usage

```
nextflow run main.nf \
--bam_dir \  path to directory of bam files
--input_suffix \  input bam file suffix
--region_name \  name or identifier for samples in batch
--region \  coordinates in chr-start-end format (Will change depending on reference genome)
--metadata_table \  list of samples
--output_dir \  publish location for outputs
--final_output_name \  final summary file name
--nested_bams \  flag that looks for bams in nested sub-directories ${bam_dir}/*/${sample_id}*${input_suffix} (default is FALSE; to run just add the flag to the nextflow command)
-resume
```

## Running on ongoing projects
For ongoing projects, run with the same output dir and final output name to have new samples append to the final summary files. 

For NSC0276_CVDCarriers, run with :
```
nextflow run main.nf \
--bam_dir /n/alignments/NSC0276_CVDCarriers \
--input_suffix 'CHM13*null-5mCG_5hmCG*.phased.bam' \
--region_name CHM13_NSC0276_CVDCarriers \
--region chrX:151389254-153479422 \
--metadata_table resources/sample_file.tsv \
--output_dir NSC0276_CVDCarriers_CHM13_samples_nested_bam_testing \
--final_output_name NSC0276_CVDCarriers \
--nested_bams \
-resume
```
For NSC0268_Maureen, run with:
```
nextflow run main.nf \
--bam_dir /n/alignments/NSC0268_Maureen \
--input_suffix 'CHM13*null-5mCG_5hmCG*.phased.bam' \
--region_name chm13_NSC0268_Maureen_samples \
--region chrX:151389254-153479422 \
--metadata_table resources/sample_file.tsv \
--output_dir NSC0268_Maureen_CHM13_samples \
--final_output_name NSC0268_Maureen_outputs \
--nested_bams \
-resume

```
### Reference genome specific coordinates:
If your starting input bam file has been aligned to the GRChg38 reference genome, then you will use chrX:153121316-155216212 as your coordinates.

If your starting input bam file has been aligned to the T2T-CHM13 reference genome, then you will use chrX:151389254-153479422 as your coordinates.

### Note about input bam directory
- If the bams are nested, use the --nested_bams flag

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
