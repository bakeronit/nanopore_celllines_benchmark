## Code and data used to generate all figures in the manuscript.

- Figure 1 was created in [bioRender](https://www.biorender.com/) under license.

- All other figures were generated using ggplot2 in R, and then assembled with either [OmniGraffle](https://www.omnigroup.com/omnigraffle) or [patchwork](https://github.com/thomasp85/patchwork).

**Figure 2**

- **R code:** `1.dna_mixing_celllines/R_analysis/1.benchmark_snv_calling.Rmd`

- **Data:** 
  - Summary of somatic SNVs and indels benchmarking results using ClairS and DeepSomatic
  - True positive and false positive calls from benchmarking.

```bash
# figure 2a
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_10.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_20.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_30.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_40.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_50.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_60.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_70.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_80.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_90.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL/summary.txt

# figure 2b
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_10.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_20.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_30.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_40.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_50.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_60.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_70.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_80.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_90.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_30.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_60.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_70.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_50.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_10.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_90.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_80.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_40.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_20.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL_passed/fp.vcf

# figure 
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_10.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_20.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_30.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_40.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_50.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_60.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_70.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_80.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_90.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL/summary.txt

#figure 2d
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_10.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_20.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_30.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_40.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_50.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_60.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_70.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_80.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_90.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_10.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_70.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_30.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_50.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_60.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_20.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_40.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_80.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_90.COLO829_BL/fp.vcf
```

**Figure 3**

- **R code:** `1.dna_mixing_celllines/R_analysis/3.benchmark_sv_calling.Rmd`

- **Data:** 
  - Merged VCF files of SV calling with short-read goldstandard for two cell lines with different purity called using nanomonsv, SAVANA, Severus, and Delly.
  
```bash
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_70.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_20.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_40.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_50.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_60.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_80.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_90.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_30.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_10.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_50.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_70.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_20.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_40.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_50.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_60.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_80.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_90.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_30.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_10.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_50.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_70.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_20.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_40.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_50.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_60.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_80.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_90.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_30.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_10.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_50.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_70.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_20.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_40.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_50.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_60.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_80.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_90.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_30.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_10.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_50.COLO829_BL.merged.vcf
```

**Figure 4**

- **R code:** `2.simulate_sequencing_depth/R_analysis/1.benchmark_snv_calling.Rmd`

- **Data:** 
  - Summary of somatic SNVs and indels benchmarking results using ClairS and DeepSomatic for samples with 100% tumour purity
  
```bash
# figure 4a
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.30x.HCC1937_BL.30x/summary.txt

# figure 4b
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.30x.HCC1937_BL.30x/summary.txt
```

**Figure 5**

- **R code:**
  - Figure 5a: `1.dna_mixing_celllines/R_analysis/1.snv_mutational_signature.Rmd`
  - Figure 5bc: `3.igv_check/R_analysis/lr_unique_sv_length_type.Rmd`

- **Data:** 
  - VCF files of false positive SNVs and short-read goldstandard SNVs.
  - BEDPE files all SVs detected in samples with 100% tumour purity labelled with suppvec indicating LR-unique, concordant, or SR-unique.
  
```bash
# figure 5a
goldstandard/vcfs/colo829/merged_normed_isec_snv.goldstandard.vcf.gz.tbi
goldstandard/vcfs/colo829/merged_normed_isec_snv.goldstandard.vcf.gz
goldstandard/vcfs/hcc1937/merged_normed_isec_snv.goldstandard.vcf.gz.tbi
goldstandard/vcfs/hcc1937/merged_normed_isec_snv.goldstandard.vcf.gz
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_30.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_60.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_70.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_50.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_10.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_90.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_80.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_40.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_20.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL_passed/fp.vcf

# figure 5bc
3.igv_check/lr_specific/bedpe/COLO829.COLO829_BL.all.bedpe
3.igv_check/lr_specific/bedpe/HCC1937.HCC1937_BL.all.bedpe
```

**Supplementary figure 1**

- **R code:**`1.dna_mixing_celllines/R_analysis/0.purity_and_qc.Rmd`

- **Data:** 
  - Data yield and read N50 for sequencing
  - Stats files for bam files generated by mosdepth, bamcov, and alignment N50

```bash
# figure s1a
1.dna_mixing_celllines/yield.txt
1.dna_mixing_celllines.N50.txt

# figure s1b
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_40.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_30.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_BL.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_40.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_10.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_80.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_50.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_80.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_90.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_70.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_20.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_50.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_10.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_90.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_60.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_30.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_BL.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_60.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_70.mosdepth.summary.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_20.mosdepth.summary.txt

1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_50.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_80.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_40.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_70.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_50.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_60.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_70.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_20.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_BL.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_BL.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_30.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_90.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_80.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_20.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_40.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_60.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_90.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_10.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_30.bamcov.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_10.bamcov.txt

1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_90.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_60.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_40.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_30.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_80.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_BL.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_10.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_20.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_80.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_50.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_40.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_30.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_50.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_70.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_70.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_60.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_20.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_90.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/COLO829_10.bamN50.txt
1.dna_mixing_celllines/work/analysis/qc/bam/R10/sup/HCC1937_BL.bamN50.txt
```

**Supplementary figure 2**

- **R code:**`1.dna_mixing_celllines/R_analysis/x.gold_standard.Rmd`

- **Data:** 
  - List of variant sites from intersection of multiple VCF files, generated by `bcftools isec`. 
  - Merged VCF files indicate the shared SV called from different tools.

```bash
goldstandard/vcfs/colo829/isec_hom100_snv_dir/sites.txt.gz
goldstandard/vcfs/colo829/isec_hom100_snv_dir/sites.txt.gz.tbi
goldstandard/vcfs/colo829/isec_indel_dir/sites.txt.gz
goldstandard/vcfs/colo829/isec_indel_dir/sites.txt.gz.tbi
goldstandard/vcfs/hcc1937/isec_hom100_snv_dir/sites.txt.gz
goldstandard/vcfs/hcc1937/isec_hom100_snv_dir/sites.txt.gz.tbi
goldstandard/vcfs/hcc1937/isec_indel_dir/sites.txt.gz
goldstandard/vcfs/hcc1937/isec_indel_dir/sites.txt.gz.tbi
goldstandard/structural_variation/analysis/svs/jasmine_merge/colo829/colo829.final_merged.vcf
goldstandard/structural_variation/analysis/svs/jasmine_merge/hcc1937/hcc1937.final_merged.vcf
```

**Supplementary figure 3**

- **R code:**
  - Supplementary figure s3a: `1.dna_mixing_celllines/R_analysis/0.purity_and_qc.Rmd`
  - Supplmentary figure s3b: `1.dna_mixing_celllines/R_analysis/2.methylation_analysis.Rmd`

- **Data:** 
  - List of SNVs in VCF files located in genomic regions with copy number=2 
  - Bed files of methylation frequency at CpG sites.
  
```bash
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937_50.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829_90.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937_40.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829_40.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937_10.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829_10.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937_20.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937_80.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937_60.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829_20.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937_90.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829_50.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937_70.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829_80.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829_60.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/HCC1937_30.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829_30.vcf
1.dna_mixing_celllines/0.purity_check/dna_mixing/vcf_files/deepsomatic/cn2/COLO829_70.vcf

1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_20.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_70.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_80.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_70.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_50.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_90.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_20.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_10.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_50.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_30.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_60.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_BL.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_30.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_40.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_60.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_BL.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/HCC1937_40.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_10.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_80.bed.gz
1.dna_mixing_celllines/work/analysis/mod/R10/sup/COLO829_90.bed.gz
```

**Supplementary figure 4**

- **R code:** `1.dna_mixing_celllines/R_analysis/1.benchmark_snv_calling.Rmd`

- **Data:** 
  - Summary of somatic SNVs benchmarking results using ClairS and DeepSomatic
  
```bash
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_10.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_70.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_30.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_50.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_60.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_20.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_40.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_80.COLO829_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL/summary.txt
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_90.COLO829_BL/summary.txt
```

**Supplementary figure 5, 6**

- **R code:** `1.dna_mixing_celllines/R_analysis/1.benchmark_snv_calling.Rmd`

- **Data:** 
  - True positive and false positive calls of samples with diffrent tumour purity from benchmarking
  
```bash
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_30.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_30.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_60.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_60.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_70.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_70.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_50.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_50.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_10.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_10.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_90.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_90.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_80.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_80.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_40.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_40.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_20.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829_20.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL_passed/tp.vcf

1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_10.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_30.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_20.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_10.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_70.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_90.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_30.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_50.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_80.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_60.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_50.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_20.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_40.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_40.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_80.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_70.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/HCC1937_60.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/deepsomatic/COLO829_90.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_10.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_10.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_30.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_20.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_10.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_70.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_70.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_90.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_30.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_30.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_50.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_50.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_80.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_60.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_60.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_50.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_20.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_20.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_40.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_40.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_40.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_80.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_80.COLO829_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_70.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/HCC1937_60.HCC1937_BL/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_90.COLO829_BL/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/indels/somatic/R10/sup/clairS/COLO829_90.COLO829_BL/tp.vcf
```

**Supplementary figure 7**

- **R code:** `1.dna_mixing_celllines/R_analysis/3.benchmark_sv_calling.Rmd`

- **Data:** 
  - Merged SV VCF with short-read gold standard of samples with different tumour purity.
  
```bash
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_70.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_20.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_40.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_50.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_60.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_80.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_90.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_30.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937_10.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829_50.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_70.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_20.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_40.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_50.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_60.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_80.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_90.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_30.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937_10.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829_50.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_70.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_20.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_40.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_50.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_60.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_80.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_90.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_30.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937_10.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829_50.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_70.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_20.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_40.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_50.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_60.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_80.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_90.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_30.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937_10.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829_50.COLO829_BL.merged.vcf
```

**Supplementary figure 8**

- **R code:** `2.simulate_sequencing_depth/R_analysis/1.benchmark_snv_calling.Rmd`

- **Data:** 
  - True positive and false positive calls for different sequencing depth combinations and 4 tumour purity levels from benchmarking.

```bash
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.30x.COLO829_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.30x.COLO829_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.45x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.45x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.45x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.45x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.60x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.60x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.60x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.60x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.45x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.45x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.30x.COLO829_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.30x.COLO829_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.15x.HCC1937_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.15x.HCC1937_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.15x.COLO829_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.15x.COLO829_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.45x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.45x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.30x.COLO829_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.30x.COLO829_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.30x.HCC1937_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.30x.HCC1937_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.45x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.45x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.45x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.45x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.60x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.60x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.60x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.60x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.45x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.45x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.30x.COLO829_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.30x.COLO829_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.15x.HCC1937_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.15x.HCC1937_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.15x.COLO829_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.15x.COLO829_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.45x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.45x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.30x.HCC1937_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.30x.HCC1937_BL.15x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.30x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.30x_passed/tp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.15x_passed/fp.vcf
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.15x_passed/tp.vcf
```

**Supplementary figure 9**

- **R code:** `2.simulate_sequencing_depth/R_analysis/1.benchmark_snv_calling.md`

- **Data:** 
  - Summary of somatic SNVs and indels benchmarking results using ClairS and DeepSomatic for 9 sequencing depth combinations in samples with 100, 80, 60 and 40% tumour purity.

```bash
# figure s9a
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_40.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_40.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_80.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_40.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_40.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_60.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_80.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_40.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_80.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_80.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_80.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_80.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_60.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_60.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_40.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_60.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_60.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_80.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_80.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_60.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_80.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_60.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_80.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_60.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_40.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_40.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_40.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_80.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_80.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_80.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_60.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_60.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_40.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_80.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_40.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_80.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_40.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_40.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_80.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_60.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_60.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_40.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_60.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_80.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_60.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_40.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_40.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_60.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_60.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_80.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_40.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_60.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/COLO829_60.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/deepsomatic/HCC1937_40.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_40.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_40.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_80.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_40.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_40.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_60.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_80.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_40.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_80.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_80.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_80.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_80.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_60.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_60.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_40.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_60.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_60.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_80.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_80.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_60.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_80.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_60.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_80.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_60.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_40.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_40.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_40.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_80.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_80.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_80.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_60.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_60.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_40.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_80.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_40.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_80.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_40.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_40.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_80.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_60.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_60.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_40.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_60.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_80.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_60.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_40.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_40.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_60.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_60.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_80.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_40.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_60.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/COLO829_60.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS/HCC1937_40.30x.HCC1937_BL.30x/summary.txt

# figure s9b
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_40.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_40.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_80.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_40.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_40.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_60.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_80.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_40.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_80.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_80.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_80.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_80.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_60.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_60.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_40.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_60.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_60.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_80.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_80.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_60.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_80.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_60.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_80.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_60.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_40.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_40.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_40.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_80.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_80.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_80.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_60.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_60.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_40.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_80.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_40.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_80.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_40.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_40.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_80.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_60.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_60.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_40.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_60.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_80.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_60.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_40.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_40.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_60.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_60.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_80.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_40.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_60.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/COLO829_60.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/deepsomatic/HCC1937_40.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_40.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_40.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_80.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_40.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_40.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_60.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_80.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_40.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_80.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_80.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_80.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_80.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_60.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_60.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_40.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_60.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_60.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_80.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_80.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_60.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_80.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_60.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_80.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_60.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.60x.HCC1937_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_40.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.15x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_40.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_40.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_80.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_80.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_80.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_60.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_60.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_40.45x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_80.60x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_40.45x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_80.60x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_40.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_40.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_80.45x.HCC1937_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_60.45x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_60.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_40.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_60.60x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.30x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_80.60x.COLO829_BL.60x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_60.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_40.15x.HCC1937_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937.30x.HCC1937_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_40.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_60.45x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_60.30x.COLO829_BL.15x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_80.45x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_40.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_60.30x.COLO829_BL.30x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/COLO829_60.60x.COLO829_BL.45x/summary.txt
2.simulate_sequencing_depth/analysis/benchmark/indels/somatic/clairS/HCC1937_40.30x.HCC1937_BL.30x/summary.txt
```

**Supplementary figure 10**

- **R code:** `2.simulate_sequencing_depth/R_analysis/2.benchmark_sv_calling.Rmd`

- **Data:** 
  - Merged VCF files of SV calling for different sequencing depth combinations in samples with 100% tumour purity.
  
```bash
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/COLO829.45x.COLO829_BL.30x/COLO829.45x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/COLO829.30x.COLO829_BL.15x/COLO829.30x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/HCC1937.15x.HCC1937_BL.15x/HCC1937.15x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/HCC1937.45x.HCC1937_BL.15x/HCC1937.45x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/COLO829.45x.COLO829_BL.45x/COLO829.45x.COLO829_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/COLO829.45x.COLO829_BL.15x/COLO829.45x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/HCC1937.60x.HCC1937_BL.30x/HCC1937.60x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/HCC1937.45x.HCC1937_BL.30x/HCC1937.45x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/HCC1937.60x.HCC1937_BL.45x/HCC1937.60x.HCC1937_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/COLO829.60x.COLO829_BL.30x/COLO829.60x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/HCC1937.60x.HCC1937_BL.60x/HCC1937.60x.HCC1937_BL.60x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/COLO829.15x.COLO829_BL.15x/COLO829.15x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/COLO829.60x.COLO829_BL.45x/COLO829.60x.COLO829_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/COLO829.60x.COLO829_BL.60x/COLO829.60x.COLO829_BL.60x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/HCC1937.45x.HCC1937_BL.45x/HCC1937.45x.HCC1937_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/HCC1937.30x.HCC1937_BL.15x/HCC1937.30x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/COLO829.30x.COLO829_BL.30x/COLO829.30x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/nanomonsv/HCC1937.30x.HCC1937_BL.30x/HCC1937.30x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/COLO829.45x.COLO829_BL.30x/COLO829.45x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/COLO829.30x.COLO829_BL.15x/COLO829.30x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/HCC1937.15x.HCC1937_BL.15x/HCC1937.15x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/HCC1937.45x.HCC1937_BL.15x/HCC1937.45x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/COLO829.45x.COLO829_BL.45x/COLO829.45x.COLO829_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/COLO829.45x.COLO829_BL.15x/COLO829.45x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/HCC1937.60x.HCC1937_BL.30x/HCC1937.60x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/HCC1937.45x.HCC1937_BL.30x/HCC1937.45x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/HCC1937.60x.HCC1937_BL.45x/HCC1937.60x.HCC1937_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/COLO829.60x.COLO829_BL.30x/COLO829.60x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/HCC1937.60x.HCC1937_BL.60x/HCC1937.60x.HCC1937_BL.60x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/COLO829.15x.COLO829_BL.15x/COLO829.15x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/COLO829.60x.COLO829_BL.45x/COLO829.60x.COLO829_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/COLO829.60x.COLO829_BL.60x/COLO829.60x.COLO829_BL.60x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/HCC1937.45x.HCC1937_BL.45x/HCC1937.45x.HCC1937_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/HCC1937.30x.HCC1937_BL.15x/HCC1937.30x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/COLO829.30x.COLO829_BL.30x/COLO829.30x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/delly/HCC1937.30x.HCC1937_BL.30x/HCC1937.30x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/COLO829.45x.COLO829_BL.30x/COLO829.45x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/COLO829.30x.COLO829_BL.15x/COLO829.30x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/HCC1937.15x.HCC1937_BL.15x/HCC1937.15x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/HCC1937.45x.HCC1937_BL.15x/HCC1937.45x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/COLO829.45x.COLO829_BL.45x/COLO829.45x.COLO829_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/COLO829.45x.COLO829_BL.15x/COLO829.45x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/HCC1937.60x.HCC1937_BL.30x/HCC1937.60x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/HCC1937.45x.HCC1937_BL.30x/HCC1937.45x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/HCC1937.60x.HCC1937_BL.45x/HCC1937.60x.HCC1937_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/COLO829.60x.COLO829_BL.30x/COLO829.60x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/HCC1937.60x.HCC1937_BL.60x/HCC1937.60x.HCC1937_BL.60x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/COLO829.15x.COLO829_BL.15x/COLO829.15x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/COLO829.60x.COLO829_BL.45x/COLO829.60x.COLO829_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/COLO829.60x.COLO829_BL.60x/COLO829.60x.COLO829_BL.60x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/HCC1937.45x.HCC1937_BL.45x/HCC1937.45x.HCC1937_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/HCC1937.30x.HCC1937_BL.15x/HCC1937.30x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/COLO829.30x.COLO829_BL.30x/COLO829.30x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/severus/HCC1937.30x.HCC1937_BL.30x/HCC1937.30x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/COLO829.45x.COLO829_BL.30x/COLO829.45x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/COLO829.30x.COLO829_BL.15x/COLO829.30x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/HCC1937.15x.HCC1937_BL.15x/HCC1937.15x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/HCC1937.45x.HCC1937_BL.15x/HCC1937.45x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/COLO829.45x.COLO829_BL.45x/COLO829.45x.COLO829_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/COLO829.45x.COLO829_BL.15x/COLO829.45x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/HCC1937.60x.HCC1937_BL.30x/HCC1937.60x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/HCC1937.45x.HCC1937_BL.30x/HCC1937.45x.HCC1937_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/HCC1937.60x.HCC1937_BL.45x/HCC1937.60x.HCC1937_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/COLO829.60x.COLO829_BL.30x/COLO829.60x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/HCC1937.60x.HCC1937_BL.60x/HCC1937.60x.HCC1937_BL.60x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/COLO829.15x.COLO829_BL.15x/COLO829.15x.COLO829_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/COLO829.60x.COLO829_BL.45x/COLO829.60x.COLO829_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/COLO829.60x.COLO829_BL.60x/COLO829.60x.COLO829_BL.60x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/HCC1937.45x.HCC1937_BL.45x/HCC1937.45x.HCC1937_BL.45x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/HCC1937.30x.HCC1937_BL.15x/HCC1937.30x.HCC1937_BL.15x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/COLO829.30x.COLO829_BL.30x/COLO829.30x.COLO829_BL.30x.merged.vcf
2.simulate_sequencing_depth/analysis/benchmark/svs/savana/HCC1937.30x.HCC1937_BL.30x/HCC1937.30x.HCC1937_BL.30x.merged.vcf
```

**Supplementary figure 11, 12, 13, 14, 15**

- **R code:** `1.dna_mixing_celllines/R_analysis/4.genome_regions.Rmd`

- **Data:** 
  - True positive and false positive calls of samples with 100% tumour purity from benchmarking

```bash
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/COLO829.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/clairS/HCC1937.HCC1937_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/COLO829.COLO829_BL_passed/tp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL_passed/fp.vcf
1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup/deepsomatic/HCC1937.HCC1937_BL_passed/tp.vcf
```

**Supplementary figure 16**

- **R code:** `1.dna_mixing_celllines/R_analysis/4.genome_regions.Rmd`

- **Data:** 
  - Merged SV VCF of samples with 100% tumour purity from comparing with short-read goldstandard.

```bash
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/nanomonsv/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/delly/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/severus/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/HCC1937.HCC1937_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark/svs/savana/R10/sup/COLO829.COLO829_BL.merged.vcf
```

**Supplementary figure 17**

- **R code:** `3.igv_check/R_analysis/germline_leakage.Rmd`

- **Data:** 
  - Tabular files with each SNV mutation annotated with gnoMAD AF_grpmax values

```bash
#figure s17a
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937_30.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937_20.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937_10.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/COLO829_70.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937_90.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/COLO829.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/COLO829_30.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/COLO829_50.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937_80.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/COLO829_60.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937_50.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/COLO829_20.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937_40.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/COLO829_40.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/COLO829_80.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937_70.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/HCC1937_60.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/deepsomatic/R10/sup/COLO829_90.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829_10.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937_30.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937_20.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937_10.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829_70.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937_90.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829_30.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829_50.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937_80.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829_60.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937_50.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829_20.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937_40.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829_40.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829_80.COLO829_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937_70.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/HCC1937_60.HCC1937_BL/snv_gnomad_af_anno.tsv
1.dna_mixing_celllines/work/analysis/snvs/clairS/R10/sup/COLO829_90.COLO829_BL/snv_gnomad_af_anno.tsv

#figure s17b
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/COLO829.45x.COLO829_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/COLO829.30x.COLO829_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/HCC1937.15x.HCC1937_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/HCC1937.45x.HCC1937_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/COLO829.45x.COLO829_BL.45x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/COLO829.45x.COLO829_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/HCC1937.60x.HCC1937_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/HCC1937.45x.HCC1937_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/HCC1937.60x.HCC1937_BL.45x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/COLO829.60x.COLO829_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/HCC1937.60x.HCC1937_BL.60x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/COLO829.15x.COLO829_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/COLO829.60x.COLO829_BL.45x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/COLO829.60x.COLO829_BL.60x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/HCC1937.45x.HCC1937_BL.45x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/HCC1937.30x.HCC1937_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/COLO829.30x.COLO829_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/deepsomatic/HCC1937.30x.HCC1937_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/COLO829.45x.COLO829_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/COLO829.30x.COLO829_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/HCC1937.15x.HCC1937_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/HCC1937.45x.HCC1937_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/COLO829.45x.COLO829_BL.45x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/COLO829.45x.COLO829_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/HCC1937.60x.HCC1937_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/HCC1937.45x.HCC1937_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/HCC1937.60x.HCC1937_BL.45x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/COLO829.60x.COLO829_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/HCC1937.60x.HCC1937_BL.60x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/COLO829.15x.COLO829_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/COLO829.60x.COLO829_BL.45x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/COLO829.60x.COLO829_BL.60x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/HCC1937.45x.HCC1937_BL.45x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/HCC1937.30x.HCC1937_BL.15x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/COLO829.30x.COLO829_BL.30x/snv_gnomad_af_anno.tsv
2.simulate_sequencing_depth/analysis/snvs/clairS/HCC1937.30x.HCC1937_BL.30x/snv_gnomad_af_anno.tsv
```

**Supplementary figure 18**

- **R code:** `1.dna_mixing_celllines/R_analysis/x.gold_standard.Rmd`

- **Data:** 
  - Merged SV VCF files with public COLO829 SV truthset.

```bash
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/nanomonsv/R10/sup/COLO829_50.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/delly/R10/sup/COLO829_50.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/severus/R10/sup/COLO829_50.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829_60.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829_70.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829_40.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829_20.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829_30.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829_10.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829_90.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829_80.COLO829_BL.merged.vcf
1.dna_mixing_celllines/work/analysis/benchmark_public/svs/savana/R10/sup/COLO829_50.COLO829_BL.merged.vcf
```
