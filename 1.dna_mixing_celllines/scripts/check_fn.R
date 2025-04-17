library(tidyverse)

df <- read_tsv("fn.vcf",col_names = c("chrom","pos","id","ref","alt","qual","filter","info","format")) 
df |> filter(chrom %in% paste0("chr",c(seq(1:22),"X","Y"))) |> 
  mutate(chrom=factor(chrom, levels=paste0("chr",c(seq(1:22),"X","Y")))) |> 
  ggplot(aes(x=chrom)) + geom_bar(stat="count") + labs(x="")


df |> select(chrom, pos,ref, alt) |> filter(nchar(alt)==1) |>  unite(snv,ref,alt,sep=">") |> 
  mutate(snv=case_when(
    snv == "A>C" ~ "T>G",
    snv == "A>T" ~ "T>A",
    snv == "A>G" ~ "T>C",
    snv == "G>A" ~ "C>T",
    snv == "G>C" ~ "C>G",
    snv == "G>T" ~ "C>A",
    .default = snv
  )) |> 
  ggplot(aes(x=snv)) + geom_bar(stat="count") + 
  geom_text(aes(label = ..count..), stat = "count",hjust=1.5,colour = "white") + 
  coord_flip() + labs(x="")


