library(RColorBrewer)
library(ggpubfigs)

workdir <- "/working/lab_nicw/jiaZ/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/"
genome_stratification_dir <- "/mnt/backedup/home/jiaZ/working/data/genome-stratification/"
gs_dir <- "/mnt/backedup/home/jiaZ/working/general/goldstandard/"

purity_pals <- c(brewer.pal(9,"YlOrRd"), "#380010")
ggsci_pals <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF")

colorblind_pals <- friendly_pals$contrast_three

public_tool_names = c("clairS" = "ClairS", "deepsomatic" = "DeepSomatic", "delly"="Delly","nanomonsv" = "nanomonsv", "savana" = "SAVANA", "severus"= "Severus")
valid_chr <- paste0("chr",c(1:22,"X","Y"))