R Session Info
================

All analyses were performed using **R 4.3.2** within an Rstudio Docker
container
([bakeronit/rstudio_hpc_cancer:0.1.3](https://hub.docker.com/layers/bakeronit/rstudio_hpc_cancer/0.1.3/images/sha256-e288941e0ab0b4a5e2152225169853912342132fbae370241e1d26912b270d5a)).
A complete list of used package versions and session details is provided
below to ensure reproducibility.

    ## R version 4.3.2 (2023-10-31)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 22.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Etc/UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] introdataviz_0.0.0.9003                   
    ##  [2] circlize_0.4.15                           
    ##  [3] VariantAnnotation_1.48.0                  
    ##  [4] Rsamtools_2.18.0                          
    ##  [5] UpSetR_1.4.0                              
    ##  [6] gtools_3.9.5                              
    ##  [7] IlluminaHumanMethylationEPICmanifest_0.3.0
    ##  [8] minfi_1.48.0                              
    ##  [9] bumphunter_1.44.0                         
    ## [10] locfit_1.5-9.8                            
    ## [11] iterators_1.0.14                          
    ## [12] foreach_1.5.2                             
    ## [13] SummarizedExperiment_1.32.0               
    ## [14] Biobase_2.62.0                            
    ## [15] MatrixGenerics_1.14.0                     
    ## [16] matrixStats_1.1.0                         
    ## [17] ggsci_3.0.0                               
    ## [18] RColorBrewer_1.1-3                        
    ## [19] deconstructSigs_1.8.0                     
    ## [20] ggpubfigs_0.0.1                           
    ## [21] ggridges_0.5.4                            
    ## [22] BSgenome.Hsapiens.UCSC.hg38_1.4.5         
    ## [23] BSgenome_1.70.1                           
    ## [24] rtracklayer_1.62.0                        
    ## [25] BiocIO_1.12.0                             
    ## [26] Biostrings_2.70.1                         
    ## [27] XVector_0.42.0                            
    ## [28] GenomicRanges_1.54.1                      
    ## [29] GenomeInfoDb_1.38.5                       
    ## [30] IRanges_2.36.0                            
    ## [31] S4Vectors_0.40.2                          
    ## [32] BiocGenerics_0.48.1                       
    ## [33] patchwork_1.1.3                           
    ## [34] ggh4x_0.2.8                               
    ## [35] data.table_1.14.10                        
    ## [36] lubridate_1.9.3                           
    ## [37] forcats_1.0.0                             
    ## [38] stringr_1.5.0                             
    ## [39] dplyr_1.1.3                               
    ## [40] purrr_1.0.2                               
    ## [41] readr_2.1.4                               
    ## [42] tidyr_1.3.0                               
    ## [43] tibble_3.2.1                              
    ## [44] ggplot2_3.4.4                             
    ## [45] tidyverse_2.0.0                           
    ## [46] here_1.0.1                                
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] shape_1.4.6               rstudioapi_0.15.0        
    ##   [3] magrittr_2.0.3            GenomicFeatures_1.54.1   
    ##   [5] rmarkdown_2.25            GlobalOptions_0.1.2      
    ##   [7] zlibbioc_1.48.0           vctrs_0.6.4              
    ##   [9] multtest_2.58.0           memoise_2.0.1            
    ##  [11] DelayedMatrixStats_1.24.0 RCurl_1.98-1.13          
    ##  [13] askpass_1.2.0             htmltools_0.5.8.1        
    ##  [15] S4Arrays_1.2.0            progress_1.2.2           
    ##  [17] curl_5.1.0                Rhdf5lib_1.24.0          
    ##  [19] SparseArray_1.2.2         rhdf5_2.46.0             
    ##  [21] nor1mix_1.3-0             plyr_1.8.9               
    ##  [23] cachem_1.0.8              GenomicAlignments_1.38.0 
    ##  [25] lifecycle_1.0.4           pkgconfig_2.0.3          
    ##  [27] Matrix_1.6-1.1            R6_2.5.1                 
    ##  [29] fastmap_1.1.1             GenomeInfoDbData_1.2.11  
    ##  [31] digest_0.6.33             siggenes_1.76.0          
    ##  [33] colorspace_2.1-0          reshape_0.8.9            
    ##  [35] AnnotationDbi_1.64.1      rprojroot_2.0.4          
    ##  [37] RSQLite_2.3.3             base64_2.0.1             
    ##  [39] filelock_1.0.2            fansi_1.0.5              
    ##  [41] timechange_0.2.0          httr_1.4.7               
    ##  [43] abind_1.4-5               compiler_4.3.2           
    ##  [45] beanplot_1.3.1            rngtools_1.5.2           
    ##  [47] bit64_4.0.5               withr_2.5.2              
    ##  [49] BiocParallel_1.36.0       DBI_1.1.3                
    ##  [51] HDF5Array_1.30.0          biomaRt_2.58.0           
    ##  [53] MASS_7.3-60               openssl_2.1.1            
    ##  [55] rappdirs_0.3.3            DelayedArray_0.28.0      
    ##  [57] rjson_0.2.21              tools_4.3.2              
    ##  [59] quadprog_1.5-8            glue_1.8.0               
    ##  [61] restfulr_0.0.15           nlme_3.1-163             
    ##  [63] rhdf5filters_1.14.1       grid_4.3.2               
    ##  [65] generics_0.1.3            gtable_0.3.4             
    ##  [67] tzdb_0.4.0                preprocessCore_1.64.0    
    ##  [69] hms_1.1.3                 xml2_1.3.5               
    ##  [71] utf8_1.2.4                pillar_1.9.0             
    ##  [73] limma_3.58.1              genefilter_1.84.0        
    ##  [75] splines_4.3.2             BiocFileCache_2.10.1     
    ##  [77] lattice_0.21-9            survival_3.5-7           
    ##  [79] bit_4.0.5                 GEOquery_2.70.0          
    ##  [81] annotate_1.80.0           tidyselect_1.2.1         
    ##  [83] knitr_1.45                gridExtra_2.3            
    ##  [85] xfun_0.41                 scrime_1.3.5             
    ##  [87] statmod_1.5.0             stringi_1.7.12           
    ##  [89] yaml_2.3.7                evaluate_0.23            
    ##  [91] codetools_0.2-19          cli_3.6.1                
    ##  [93] xtable_1.8-4              munsell_0.5.0            
    ##  [95] Rcpp_1.0.12               dbplyr_2.4.0             
    ##  [97] png_0.1-8                 XML_3.99-0.15            
    ##  [99] blob_1.2.4                prettyunits_1.2.0        
    ## [101] mclust_6.0.0              doRNG_1.8.6              
    ## [103] sparseMatrixStats_1.14.0  bitops_1.0-7             
    ## [105] illuminaio_0.44.0         scales_1.3.0             
    ## [107] crayon_1.5.2              rlang_1.1.5              
    ## [109] KEGGREST_1.42.0
