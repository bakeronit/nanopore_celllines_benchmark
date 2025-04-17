library(tidyverse)
library(RColorBrewer)
mypal <- c(brewer.pal(9,"YlOrRd"),"black")

# check the in silico mixing

df1 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829.vcf", col_names = F,comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="100")
info_cols <- df1$X9 |> str_split(":") |> first()
df1 <- df1 |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="100")
df2 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829-0.9.vcf",col_names = F, comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="90")
df3 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829-0.8.vcf",col_names = F, comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="80")
df4 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829-0.7.vcf",col_names = F, comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="70")
df5 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829-0.6.vcf",col_names = F, comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="60")
df6 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829-0.5.vcf",col_names = F, comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="50")
df7 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829-0.4.vcf",col_names = F, comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="40")
df8 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829-0.3.vcf",col_names = F, comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="30")
df9 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829-0.2.vcf",col_names = F, comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="20")
df10 <- read_tsv("working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/0.purity_check/in_silico/vcf_files/cn2/COLO829-0.1.vcf",col_names = F, comment = "#") |> separate_wider_delim(cols = "X10",delim = ":",names = info_cols) |> add_column(purity="10")

rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10) |> mutate(purity=factor(purity,levels=seq(10,100,10) |> as.character())) |> 
  ggplot(aes(x=as.double(VAF),color=purity)) + geom_density(linewidth=1) + 
  scale_color_manual(values=mypal) + theme_bw(base_size = 12) + 
  labs(x="VAF",y="Density",color="Purity",title = "COLO829 mixed in silico") + 
  theme(legend.position = "none")

