
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

- [Base calling and read alignment for cell line mixtures](README.md)
- [Sequencing depth combinations](README.md)
- [Variant calling](README.md)

### Downstream analysis and plots

- [Construct gold
  standard](1.dna_mixing_cellines/R_analysis/x.gold_standard.md): Rmd
  file [gold_standard.Rmd](1.dna_mixing_celllines/x.gold_standard.Rmd)
- [Tumour purity affects SNV and indel
  calling](1.dna_mixing_cellines/R_analysis/1.benchmark_snv_calling.md):
  Rmd file
  [benchmark_snv_calling.Rmd](1.dna_mixing_cellines/R_analysis/1.benchmark_snv_calling.Rmd)
- [Tumour purity affects SV
  calling](1.dna_mixing_cellines/R_analysis/3.benchmark_sv_calling.md):
  Rmd file
  [benchmark_sv_calling.Rmd](1.dna_mixing_cellines/R_analysis/3.benchmark_sv_calling.md)
- [Sequencing depth affects SNV and indel calling]()
- [Sequencing depth affects SV calling]()
- [Genomic regions of false positive calls]()
- [Mutational signature analysis]()
- [Germline leakage]()
- [SV type and length in LR]()

### Miscellaneous

- [Sequencing depth check]()
- [Methylation analysis]()
- [Cicos plot]()

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE)
file for details.
