## Base calling and read alignment

Raw POD5 files from Nanopore Promethion sequencing was obtained from QIMR Berghofer. A complete list of all these files is provided as [sample_list.csv](config/sample_list.csv), the POD5 file is hug (66Tb) for any platform to host thus request by email the author.

Base calling was performed using [Dorado](https://github.com/nanoporetech/dorado) version [0.5.1](https://github.com/nanoporetech/dorado/releases/tag/v0.5.1) with canonical model [dna_r10.4.1_e8.2_400bps_sup@v4.3.0](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models) and remora model [dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mCG_5hmCG@v1](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models). The command line used was described in [dorado_basecalling.smk](workflow/rules/dorado_basecalling.smk).


Raw sequence data was stored as unaligned BAM files for each flowcell, reads with `qs<10` were excluded and converted into FASTQ format using SAMtools then aligned with minimap2, alignments of flowcells for each sample were merged into one BAM file.. The Snakemake workflow file provides details of steps and parameters used: [minimap2_align.smk](workflow/rules/minimap2_align.smk).


## Variant calling

Somatic SNV, indels and SV were called for 20 tumour samples with 10-100% tumour purity of two cell lines using the data from correspond B lymphoblast (BL) cell line as normal sample. For read depth simulation, BL sample with different read depth was used normal results in 72 combinations in each variant calling. 

BAM files for all samples were indexed with `samtools`. For soamtic SNV & indel calling, DeepSomatic (v1.6.0) and ClairS () were run as in [deepsomatic.smk](workflow/rules/snv_calling/deepsomatic.smk) and  [clairS.smk](workflow/rules/snv_calling/clairS.smk).

For benchmarking, the raw variants with PASS filter were used.