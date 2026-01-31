Supplementary note: simulation with nanoSim
================

We implemented a small in-silico simulation with NanoSim from synthetic
tumour genome, to provide a controlled stress-test of somatic variant
calling (SNVs) across tumour purity and sequencing depth while keeping
the main benchmarking focused on empirical ONT mixture datasets.

# Overview

## Goal

This simulation is designed as a *supplementary sanity check*
(orthogonal to the real mixture benchmarking) to isolate expected
failure modes as a function of:

- Tumour purity (fraction of tumour reads in the mixture)
- Sequencing depth (15×–60× on chr22)

We simulate **matched tumour/normal** ONT datasets on a **single
chromosome (chr22)** with: - A shared germline background present in
both tumour and normal - Additional somatic SNVs present only in tumour
(the ground truth)

## Inputs

- `chr22.fasta` from GRCh38 (reference genome FASTA; indexed)
- A real ONT FASTQ for NanoSim training, COLO829_BL.bam
- Synthetic germline genome and tumour genome for COLO829_BL and
  COLO829.

**Tools:**

- `SAMtools`
- `BCFtools`
- `NanoSim`
- `cstag`
- `chopper`
- `minimap2`
- `clairS`
- `DeepSomatic`

See detailed steps in the [Snakefile](../Snakefile)

## Bechmarking results

| Tool        | Sample      | Tumour purity | Tumour depth | Normal depth | Precision | Recall | F1.score |
|:------------|:------------|--------------:|:-------------|:-------------|----------:|-------:|---------:|
| clairS      | COLO829_100 |           100 | 60x          | 60x          |    0.9982 | 1.0000 |   0.9991 |
| clairS      | COLO829_100 |           100 | 60x          | 45x          |    0.9982 | 1.0000 |   0.9991 |
| clairS      | COLO829_100 |           100 | 60x          | 30x          |    0.9982 | 1.0000 |   0.9991 |
| clairS      | COLO829_100 |           100 | 45x          | 45x          |    0.9448 | 1.0000 |   0.9716 |
| clairS      | COLO829_100 |           100 | 45x          | 30x          |    0.9320 | 1.0000 |   0.9648 |
| clairS      | COLO829_100 |           100 | 45x          | 15x          |    0.9432 | 1.0000 |   0.9708 |
| clairS      | COLO829_100 |           100 | 30x          | 30x          |    0.8562 | 1.0000 |   0.9226 |
| clairS      | COLO829_100 |           100 | 30x          | 15x          |    0.8810 | 1.0000 |   0.9368 |
| clairS      | COLO829_100 |           100 | 15x          | 15x          |    0.9525 | 0.9872 |   0.9695 |
| clairS      | COLO829_80  |            80 | 60x          | 60x          |    1.0000 | 1.0000 |   1.0000 |
| clairS      | COLO829_80  |            80 | 60x          | 45x          |    0.9982 | 1.0000 |   0.9991 |
| clairS      | COLO829_80  |            80 | 60x          | 30x          |    0.9964 | 1.0000 |   0.9982 |
| clairS      | COLO829_80  |            80 | 45x          | 45x          |    0.9734 | 1.0000 |   0.9865 |
| clairS      | COLO829_80  |            80 | 45x          | 30x          |    0.9288 | 1.0000 |   0.9631 |
| clairS      | COLO829_80  |            80 | 45x          | 15x          |    0.9383 | 0.9982 |   0.9673 |
| clairS      | COLO829_80  |            80 | 30x          | 30x          |    0.9058 | 1.0000 |   0.9506 |
| clairS      | COLO829_80  |            80 | 30x          | 15x          |    0.8722 | 0.9964 |   0.9302 |
| clairS      | COLO829_80  |            80 | 15x          | 15x          |    0.9749 | 0.9927 |   0.9837 |
| clairS      | COLO829_60  |            60 | 60x          | 60x          |    1.0000 | 1.0000 |   1.0000 |
| clairS      | COLO829_60  |            60 | 60x          | 45x          |    0.9964 | 1.0000 |   0.9982 |
| clairS      | COLO829_60  |            60 | 60x          | 30x          |    0.9964 | 1.0000 |   0.9982 |
| clairS      | COLO829_60  |            60 | 45x          | 45x          |    0.9892 | 1.0000 |   0.9946 |
| clairS      | COLO829_60  |            60 | 45x          | 30x          |    0.9320 | 1.0000 |   0.9648 |
| clairS      | COLO829_60  |            60 | 45x          | 15x          |    0.9545 | 0.9964 |   0.9750 |
| clairS      | COLO829_60  |            60 | 30x          | 30x          |    0.9682 | 1.0000 |   0.9838 |
| clairS      | COLO829_60  |            60 | 30x          | 15x          |    0.8962 | 0.9927 |   0.9420 |
| clairS      | COLO829_60  |            60 | 15x          | 15x          |    0.9852 | 0.9745 |   0.9798 |
| clairS      | COLO829_40  |            40 | 60x          | 60x          |    1.0000 | 0.9982 |   0.9991 |
| clairS      | COLO829_40  |            40 | 60x          | 45x          |    1.0000 | 0.9982 |   0.9991 |
| clairS      | COLO829_40  |            40 | 60x          | 30x          |    0.9928 | 1.0000 |   0.9964 |
| clairS      | COLO829_40  |            40 | 45x          | 45x          |    0.9945 | 0.9982 |   0.9964 |
| clairS      | COLO829_40  |            40 | 45x          | 30x          |    0.9399 | 0.9982 |   0.9681 |
| clairS      | COLO829_40  |            40 | 45x          | 15x          |    0.9199 | 0.9854 |   0.9515 |
| clairS      | COLO829_40  |            40 | 30x          | 30x          |    0.9891 | 0.9945 |   0.9918 |
| clairS      | COLO829_40  |            40 | 30x          | 15x          |    0.8805 | 0.9818 |   0.9284 |
| clairS      | COLO829_40  |            40 | 15x          | 15x          |    0.9958 | 0.8577 |   0.9216 |
| clairS      | COLO829_20  |            20 | 60x          | 60x          |    1.0000 | 0.7464 |   0.8548 |
| clairS      | COLO829_20  |            20 | 60x          | 45x          |    0.9931 | 0.7865 |   0.8778 |
| clairS      | COLO829_20  |            20 | 60x          | 30x          |    0.9893 | 0.8412 |   0.9093 |
| clairS      | COLO829_20  |            20 | 45x          | 45x          |    1.0000 | 0.8011 |   0.8896 |
| clairS      | COLO829_20  |            20 | 45x          | 30x          |    0.9257 | 0.8412 |   0.8815 |
| clairS      | COLO829_20  |            20 | 45x          | 15x          |    0.9204 | 0.8230 |   0.8690 |
| clairS      | COLO829_20  |            20 | 30x          | 30x          |    0.9977 | 0.7974 |   0.8864 |
| clairS      | COLO829_20  |            20 | 30x          | 15x          |    0.8758 | 0.7591 |   0.8133 |
| clairS      | COLO829_20  |            20 | 15x          | 15x          |    1.0000 | 0.4106 |   0.5821 |
| clairS      | COLO829_10  |            10 | 60x          | 60x          |    1.0000 | 0.1697 |   0.2902 |
| clairS      | COLO829_10  |            10 | 60x          | 45x          |    0.9832 | 0.2135 |   0.3508 |
| clairS      | COLO829_10  |            10 | 60x          | 30x          |    0.9663 | 0.3139 |   0.4738 |
| clairS      | COLO829_10  |            10 | 45x          | 45x          |    1.0000 | 0.3102 |   0.4735 |
| clairS      | COLO829_10  |            10 | 45x          | 30x          |    0.8216 | 0.4033 |   0.5410 |
| clairS      | COLO829_10  |            10 | 45x          | 15x          |    0.8516 | 0.3978 |   0.5423 |
| clairS      | COLO829_10  |            10 | 30x          | 30x          |    1.0000 | 0.3522 |   0.5209 |
| clairS      | COLO829_10  |            10 | 30x          | 15x          |    0.7262 | 0.3339 |   0.4575 |
| clairS      | COLO829_10  |            10 | 15x          | 15x          |    1.0000 | 0.1332 |   0.2351 |
| deepsomatic | COLO829_100 |           100 | 60x          | 60x          |    0.9928 | 1.0000 |   0.9964 |
| deepsomatic | COLO829_100 |           100 | 60x          | 45x          |    1.0000 | 0.9982 |   0.9991 |
| deepsomatic | COLO829_100 |           100 | 60x          | 30x          |    1.0000 | 0.9964 |   0.9982 |
| deepsomatic | COLO829_100 |           100 | 45x          | 45x          |    1.0000 | 0.9982 |   0.9991 |
| deepsomatic | COLO829_100 |           100 | 45x          | 30x          |    1.0000 | 0.9945 |   0.9973 |
| deepsomatic | COLO829_100 |           100 | 45x          | 15x          |    1.0000 | 0.9635 |   0.9814 |
| deepsomatic | COLO829_100 |           100 | 30x          | 30x          |    1.0000 | 0.9927 |   0.9963 |
| deepsomatic | COLO829_100 |           100 | 30x          | 15x          |    1.0000 | 0.9635 |   0.9814 |
| deepsomatic | COLO829_100 |           100 | 15x          | 15x          |    1.0000 | 0.9635 |   0.9814 |
| deepsomatic | COLO829_80  |            80 | 60x          | 60x          |    1.0000 | 1.0000 |   1.0000 |
| deepsomatic | COLO829_80  |            80 | 60x          | 45x          |    1.0000 | 1.0000 |   1.0000 |
| deepsomatic | COLO829_80  |            80 | 60x          | 30x          |    1.0000 | 0.9982 |   0.9991 |
| deepsomatic | COLO829_80  |            80 | 45x          | 45x          |    1.0000 | 0.9964 |   0.9982 |
| deepsomatic | COLO829_80  |            80 | 45x          | 30x          |    1.0000 | 0.9982 |   0.9991 |
| deepsomatic | COLO829_80  |            80 | 45x          | 15x          |    1.0000 | 0.9745 |   0.9871 |
| deepsomatic | COLO829_80  |            80 | 30x          | 30x          |    1.0000 | 0.9927 |   0.9963 |
| deepsomatic | COLO829_80  |            80 | 30x          | 15x          |    1.0000 | 0.9745 |   0.9871 |
| deepsomatic | COLO829_80  |            80 | 15x          | 15x          |    0.9981 | 0.9726 |   0.9852 |
| deepsomatic | COLO829_60  |            60 | 60x          | 60x          |    0.9982 | 1.0000 |   0.9991 |
| deepsomatic | COLO829_60  |            60 | 60x          | 45x          |    0.9964 | 1.0000 |   0.9982 |
| deepsomatic | COLO829_60  |            60 | 60x          | 30x          |    1.0000 | 0.9982 |   0.9991 |
| deepsomatic | COLO829_60  |            60 | 45x          | 45x          |    1.0000 | 0.9982 |   0.9991 |
| deepsomatic | COLO829_60  |            60 | 45x          | 30x          |    1.0000 | 0.9964 |   0.9982 |
| deepsomatic | COLO829_60  |            60 | 45x          | 15x          |    1.0000 | 0.9745 |   0.9871 |
| deepsomatic | COLO829_60  |            60 | 30x          | 30x          |    1.0000 | 0.9945 |   0.9973 |
| deepsomatic | COLO829_60  |            60 | 30x          | 15x          |    1.0000 | 0.9745 |   0.9871 |
| deepsomatic | COLO829_60  |            60 | 15x          | 15x          |    1.0000 | 0.9343 |   0.9660 |
| deepsomatic | COLO829_40  |            40 | 60x          | 60x          |    1.0000 | 1.0000 |   1.0000 |
| deepsomatic | COLO829_40  |            40 | 60x          | 45x          |    1.0000 | 0.9982 |   0.9991 |
| deepsomatic | COLO829_40  |            40 | 60x          | 30x          |    1.0000 | 0.9982 |   0.9991 |
| deepsomatic | COLO829_40  |            40 | 45x          | 45x          |    1.0000 | 0.9945 |   0.9973 |
| deepsomatic | COLO829_40  |            40 | 45x          | 30x          |    1.0000 | 0.9927 |   0.9963 |
| deepsomatic | COLO829_40  |            40 | 45x          | 15x          |    1.0000 | 0.9763 |   0.9880 |
| deepsomatic | COLO829_40  |            40 | 30x          | 30x          |    1.0000 | 0.9708 |   0.9852 |
| deepsomatic | COLO829_40  |            40 | 30x          | 15x          |    1.0000 | 0.9526 |   0.9757 |
| deepsomatic | COLO829_40  |            40 | 15x          | 15x          |    1.0000 | 0.7281 |   0.8427 |
| deepsomatic | COLO829_20  |            20 | 60x          | 60x          |    1.0000 | 0.9635 |   0.9814 |
| deepsomatic | COLO829_20  |            20 | 60x          | 45x          |    1.0000 | 0.9325 |   0.9651 |
| deepsomatic | COLO829_20  |            20 | 60x          | 30x          |    1.0000 | 0.9215 |   0.9592 |
| deepsomatic | COLO829_20  |            20 | 45x          | 45x          |    1.0000 | 0.8175 |   0.8996 |
| deepsomatic | COLO829_20  |            20 | 45x          | 30x          |    1.0000 | 0.7737 |   0.8724 |
| deepsomatic | COLO829_20  |            20 | 45x          | 15x          |    1.0000 | 0.7755 |   0.8736 |
| deepsomatic | COLO829_20  |            20 | 30x          | 30x          |    1.0000 | 0.5474 |   0.7075 |
| deepsomatic | COLO829_20  |            20 | 30x          | 15x          |    1.0000 | 0.5657 |   0.7226 |
| deepsomatic | COLO829_20  |            20 | 15x          | 15x          |    1.0000 | 0.2099 |   0.3469 |
| deepsomatic | COLO829_10  |            10 | 60x          | 60x          |    1.0000 | 0.4891 |   0.6569 |
| deepsomatic | COLO829_10  |            10 | 60x          | 45x          |    1.0000 | 0.4270 |   0.5985 |
| deepsomatic | COLO829_10  |            10 | 60x          | 30x          |    1.0000 | 0.3631 |   0.5328 |
| deepsomatic | COLO829_10  |            10 | 45x          | 45x          |    1.0000 | 0.2427 |   0.3906 |
| deepsomatic | COLO829_10  |            10 | 45x          | 30x          |    1.0000 | 0.2245 |   0.3666 |
| deepsomatic | COLO829_10  |            10 | 45x          | 15x          |    1.0000 | 0.2281 |   0.3715 |
| deepsomatic | COLO829_10  |            10 | 30x          | 30x          |    1.0000 | 0.0803 |   0.1486 |
| deepsomatic | COLO829_10  |            10 | 30x          | 15x          |    1.0000 | 0.1150 |   0.2062 |
| deepsomatic | COLO829_10  |            10 | 15x          | 15x          |    1.0000 | 0.0365 |   0.0704 |
