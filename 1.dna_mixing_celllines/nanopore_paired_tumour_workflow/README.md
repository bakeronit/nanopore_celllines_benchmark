## Base calling and read alignment

Raw POD5 files from Nanopore Promethion sequencing was obtained from QIMR Berghofer. A complete list of all these files is provided as [sample_list.csv](config/sample_list.csv), the POD5 file is hug (66Tb) for any platform to host thus request by email the author.

Base calling was performed using [Dorado](https://github.com/nanoporetech/dorado) version [0.5.1](https://github.com/nanoporetech/dorado/releases/tag/v0.5.1) with canonical model [dna_r10.4.1_e8.2_400bps_sup@v4.3.0](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models) and remora model [dna_r10.4.1_e8.2_400bps_sup@v4.3.0_5mCG_5hmCG@v1](https://github.com/nanoporetech/dorado?tab=readme-ov-file#dna-models). The command line used was described in [dorado_basecalling.smk](1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/dorado_basecalling.smk).


Raw sequence data was stored as unaligned BAM files for each flowcell, reads with `qs<10` were excluded and converted into FASTQ format using SAMtools then aligned with minimap2, alignments of flowcells for each sample were merged into one BAM file.. The Snakemake workflow file provides details of steps and parameters used: [minimap2_align.smk](1.dna_mixing_celllines/nanopore_paired_tumour_workflow/workflow/rules/minimap2_align.smk).