Overlapped SNV calls
================

This notebook contains analysis performed to address Reviewer \#1’s
comment:

> Weakness 1. The lack of analysis for overlapping true positives or
> false positives across methods was a missed opportunity to understand
> the potential impact of ensemble methods. Albeit, an analysis of
> overlapping calls is slightly out of scope of this work. Therefore, I
> do not see overlap analysis as a necessary revision.

All scripts we ran to generate the summary of overlapped SNVs from
ClairS and DeepSomatic are included in
[intersect_snvs.smk](../intersect_snvs.smk).

We found that large proportion of both methods are overlapped, especally
for DeepSomatic, there were not many DeepSomatic unique SNVs.

<img src="ensembl_snv_files/figure-gfm/snv_overlap-1.png" width="652.8" style="display: block; margin: auto;" />

As expected, the recall is slightly lower compared with each single
method but higher precision. In addition, the union calls have higher
recalls. In the future, using union calls and integrated with filters
which take context (allele frequency, read depth, quality scores and
even genomic regions) might be a good strategy be produce high
confidence variant calling.

| Cell line | Tumour Purity | Precision of Intersection | Recall of Intersection | Recall of Union |
|:----------|:--------------|--------------------------:|-----------------------:|----------------:|
| COLO829   | 100           |                    0.9683 |                 0.9530 |          0.9725 |
| COLO829   | 90            |                    0.9732 |                 0.9525 |          0.9722 |
| COLO829   | 80            |                    0.9750 |                 0.9457 |          0.9676 |
| COLO829   | 70            |                    0.9776 |                 0.9372 |          0.9625 |
| COLO829   | 60            |                    0.9821 |                 0.9264 |          0.9575 |
| COLO829   | 50            |                    0.9842 |                 0.9035 |          0.9451 |
| COLO829   | 40            |                    0.9873 |                 0.8673 |          0.9289 |
| COLO829   | 30            |                    0.9889 |                 0.7861 |          0.8947 |
| COLO829   | 20            |                    0.9902 |                 0.5759 |          0.7714 |
| COLO829   | 10            |                    0.9912 |                 0.1797 |          0.3661 |
| HCC1937   | 100           |                    0.8751 |                 0.9517 |          0.9756 |
| HCC1937   | 90            |                    0.8852 |                 0.9375 |          0.9655 |
| HCC1937   | 80            |                    0.8926 |                 0.9295 |          0.9643 |
| HCC1937   | 70            |                    0.8987 |                 0.9146 |          0.9582 |
| HCC1937   | 60            |                    0.9046 |                 0.8878 |          0.9412 |
| HCC1937   | 50            |                    0.9110 |                 0.8522 |          0.9223 |
| HCC1937   | 40            |                    0.9196 |                 0.7454 |          0.8746 |
| HCC1937   | 30            |                    0.9229 |                 0.6438 |          0.8317 |
| HCC1937   | 20            |                    0.9289 |                 0.4070 |          0.6861 |
| HCC1937   | 10            |                    0.9422 |                 0.1150 |          0.3400 |
