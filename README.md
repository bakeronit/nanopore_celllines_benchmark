
# Cancer genome standards for long-read sequencing using cancer cell line mixtures

## Summary

This repository hosts the analysis scripts and pipeline associated with
the paper in submission.

In this study, we evaluated the performance of long-read sequencing
(LRS) for detecting somatic variants across a range of tumor purities
and sequencing depths, comparing results to short-read sequencing. We
generated 22 whole-genome sequencing datasets from controlled mixtures
of cancer and matched normal cell lines (0%–100% tumor purity). This
design enabled benchmarking of LRS-based somatic variant detection under
realistic scenarios.

## Table of Contents

### Raw data processing

- [Base calling and read alignment for cell line
  mixtures](1.dna_mixing_celllines/nanopore_paired_tumour_workflow/README.md##Base%20calling%20and%20read%20alignment)
- [Sequencing depth combinations](2.simulate_sequencing_depth/README.md)
- [Variant
  calling](1.dna_mixing_celllines/nanopore_paired_tumour_workflow/README.md##Variant%20calling)
- [QC and purity
  check](1.dna_mixing_celllines/R_analysis/0.purity_and_qc.md): Rmd file
  [purity_and_qc.Rmd](1.dna_mixing_celllines/R_analysis/0.purity_and_qc.Rmd)

### Downstream analysis and plots

- [Construct gold
  standard](1.dna_mixing_celllines/R_analysis/x.gold_standard.md): Rmd
  file
  [gold_standard.Rmd](1.dna_mixing_celllines/R_analysis/x.gold_standard.Rmd)
- [Tumour purity affects SNV and indel
  calling](1.dna_mixing_celllines/R_analysis/1.benchmark_snv_calling.md):
  Rmd file
  [benchmark_snv_calling.Rmd](1.dna_mixing_celllines/R_analysis/1.benchmark_snv_calling.Rmd)
- [Tumour purity affects SV
  calling](1.dna_mixing_celllines/R_analysis/3.benchmark_sv_calling.md):
  Rmd file
  [benchmark_sv_calling.Rmd](1.dna_mixing_celllines/R_analysis/3.benchmark_sv_calling.Rmd)
- [Sequencing depth affects SNV and indel
  calling](2.simulate_sequencing_depth/R_analysis/1.benchmark_snv_calling.md):
  Rmd file
  [benchmark_snv_calling.Rmd](2.simulate_sequencing_depth/R_analysis/1.benchmark_snv_calling.Rmd)
- [Sequencing depth affects SV
  calling](2.simulate_sequencing_depth/R_analysis/2.benchmark_sv_calling.md):
  Rmd file
  [benchmark_sv_calling.Rmd](2.simulate_sequencing_depth/R_analysis/2.benchmark_sv_calling.Rmd)
- [Mutational signature
  analysis](1.dna_mixing_celllines/R_analysis/1.snv_mutational_signature.md):
  Rmd file
  [snv_mutational_signature.Rmd](1.dna_mixing_celllines/R_analysis/1.snv_mutational_signature.Rmd)
- [Genomic regions of
  variants](1.dna_mixing_celllines/R_analysis/4.genome_regions.md): Rmd
  file
  [genome_regions.Rmd](1.dna_mixing_celllines/R_analysis/4.genome_regions.Rmd)
- [Germline leakage against tumour purity and read
  depth](3.igv_check/R_analysis/germline_leakage.md): Rmd file
  [germline_leakage.Rmd](3.igv_check/R_analysis/germline_leakage.Rmd)
- [SV type and length in
  LRS](3.igv_check/R_analysis/lr_unique_sv_length_type.md): Rmd file
  [lr_unique_length_type.Rmd](3.igv_check/R_analysis/lr_unique_sv_length_type.Rmd)

### Miscellaneous

- [Sequencing depth
  check](2.simulate_sequencing_depth/R_analysis/0.depth_check.md): Rmd
  file
  [depth_check.Rmd](2.simulate_sequencing_depth/R_analysis/0.depth_check.Rmd)
- [Methylation
  analysis](1.dna_mixing_celllines/R_analysis/2.methylation_analysis.md):
  Rmd file
  [methylation_analysis.Rmd](1.dna_mixing_celllines/R_analysis/2.methylation_analysis.Rmd)
- [Circos plot](3.igv_check/R_analysis/0.circos_plots.md): Rmd file
  [circos_plots.Rmd](3.igv_check/R_analysis/0.circos_plots.Rmd)
- [IGV check](3.igv_check/README.md)

## Reproducibility

Each section above is available as a processed Markdown (`.md`) file.
Clicking on the links will open web-readable pages that include
explanatory text, selected commands, plots, and tables. The underlying
code used to generate these outputs is provided in the corresponding R
Markdown (`.Rmd`) files. For R environment details, see the [R Session
Info](rsessioninfo.md) page.

All necessary data required to run these notebooks will be made
available via OwnCloud. This dataset is dedicated to the public domain
under the [Creative Commons CC0 1.0 Universal (CC0 1.0) Public Domain
Dedication](https://creativecommons.org/publicdomain/zero/1.0/),
allowing unrestricted reuse.

## License

This project is licensed under the BSD 3-Clause License. See the
[LICENSE](LICENSE) file for details.
