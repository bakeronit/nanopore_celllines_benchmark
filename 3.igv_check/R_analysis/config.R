library(RColorBrewer)
library(ggpubfigs)

purity_workdir <- "1.dna_mixing_celllines/work/"
depth_workdir <- "2.simulate_sequencing_depth/"
genome_stratification_dir <- "data/genome-stratification/"
gs_dir <- "gs/goldstandard/"

purity_pals <- c(brewer.pal(9,"YlOrRd"), "#380010")

colorblind_pals <- friendly_pals$contrast_three

public_tool_names = c("clairS" = "ClairS", "deepsomatic" = "DeepSomatic", "delly"="Delly","nanomonsv" = "nanomonsv", "savana" = "SAVANA", "severus"= "Severus")
valid_chr <- paste0("chr",c(1:22,"X","Y"))
