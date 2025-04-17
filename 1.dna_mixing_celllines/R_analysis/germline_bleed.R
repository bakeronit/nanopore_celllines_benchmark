library(tidyverse)
cols<- c("chrom","pos","ref","alt","qual","vaf","depth","freq")

df <- list.files(path = "~/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/snvs", full.names = TRUE, recursive = TRUE, pattern = "^snv_gnomad_af_anno.tsv$" ) |> 
  map(\(x) read_tsv(x, col_names=cols) |> add_column(prefix=x)) |> list_rbind()

df2 <- df |> mutate(sample=str_split_i(prefix, "/",i=16), tool=str_split_i(prefix,"/",i=13)) |> 
  separate_wider_delim(sample,delim = ".", names = c("tumour","normal")) |> select(-prefix) |> 
  filter(!(chrom=="chrX" & normal=="HCC1937_BL")) |> 
  group_by(tumour,tool) |> summarise(n=n(), prop=mean(freq > 0.1, na.rm = TRUE)) |> ungroup() |> 
  mutate(celltype=ifelse(startsWith(tumour,"COLO829"),"COLO829","HCC1937"), purity=str_split_i(tumour,"_",i=2), purity=ifelse(is.na(purity),"100",purity), purity=factor(purity,levels=seq(100,10,-10)))
  


df2 |> ggplot(aes(x=purity,y=prop, color=celltype,group=celltype)) + geom_point() + geom_line() +
  facet_wrap(~tool, scale="free_y") 



depth_df <- list.files(path = "~/working/bioprojects/nanopore_celllines_benchmark/2.simulate_sequencing_depth/analysis/snvs", full.names = TRUE, recursive = TRUE, pattern = "^snv_gnomad_af_anno.tsv$" ) |> 
  map(\(x) read_tsv(x, col_names=cols) |> add_column(prefix=x)) |> list_rbind()



depth_df2 <- depth_df |> mutate(sample=str_split_i(prefix, "/",i=13), tool=str_split_i(prefix,"/",i=12)) |> 
  separate_wider_delim(sample,delim = ".", names = c("tumour","t_depth","normal","n_depth")) |> select(-prefix) |> 
  filter(!(chrom=="chrX" & normal=="HCC1937_BL")) |> 
  group_by(tumour,t_depth, n_depth, tool,) |> summarise(n=n(), prop=mean(freq > 0.1, na.rm = TRUE)) |> ungroup() |> 
  mutate(celltype=ifelse(startsWith(tumour,"COLO829"),"COLO829","HCC1937"), purity=str_split_i(tumour,"_",i=2), purity=ifelse(is.na(purity),"100",purity), purity=factor(purity,levels=seq(100,10,-10)))

depth_df2 |> 
  unite(comb,c(t_depth,n_depth),sep="-") |> 
  ggplot(aes(x=comb,y=prop,group=celltype, color=celltype)) +
  geom_point() + geom_line() + facet_wrap(purity~tool, scales = "free_y")
