library(tidyverse)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

mypal <- c(brewer.pal(9,"YlOrRd"),"black")
theme_set(theme_gray(base_size = 14))
somatic_path <- "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/1.dna_mixing_celllines/work/analysis/benchmark/snvs/somatic/R10/sup"
somatic_path_simulated <- "/mnt/backedup/home/jiaZ/working/bioprojects/nanopore_celllines_benchmark/0.pilots/mixing_experiment/in_silico/benchmarks/snvs/somatic"
summary_files <- list.files(path = somatic_path,pattern = "summary.txt",recursive = TRUE, full.names = TRUE)
summary_files_simulated <- list.files(path = somatic_path_simulated,pattern = "summary.txt",recursive = TRUE, full.names = TRUE)


df <- summary_files |> set_names(\(x) str_remove_all(x,pattern = c(paste0(somatic_path,"/"),"/summary.txt") |> paste(collapse = "|"))) |> 
  map(\(x) read.table(x,header=TRUE)) |> list_rbind(names_to = "prefix") |> separate_wider_delim(prefix, delim = "/",names=c("tool","sample")) |> 
  mutate(sample=str_split(sample,"\\.") |> map_chr(first), purity=str_split(sample,"_") |> map_chr(last)) |> 
  mutate(purity=ifelse(is.na(as.numeric(purity)), 100, as.numeric(purity))) |> 
  mutate(purity=factor(purity, levels=seq(10,100,10)))


## df for simulated data, named differently
df_simulated <-summary_files_simulated |> set_names(\(x) str_remove_all(x,pattern = c(paste0(somatic_path,"/"),"/summary.txt") |> paste(collapse = "|"))) |> 
  map(\(x) read.table(x,header=TRUE)) |> list_rbind(names_to = "prefix") |> separate_wider_delim(prefix, delim = "/",names=c("tool","sample")) |> 
  mutate( purity=str_split(sample,"-") |> map_chr(last),sample=str_split(sample,"-") |> map_chr(first), purity=as.double(purity)) |> 
  mutate(purity=purity*100) |> 
  mutate(purity=factor(purity, levels=seq(10,100,10)))
  
plot_roclike <- df |> filter(Filter=="PASS") |> 
  ggplot(aes(x=Recall,y=Precision)) + geom_path(aes(linetype=tool),linewidth=1) + 
  geom_point(size=4,aes(color=purity,shape=tool)) + scale_color_manual(values=mypal)
  
plot2 <- df |> filter(Filter=="PASS") |> pivot_longer(c(Recall,Precision)) |> 
  ggplot(aes(x=purity,y=value,fill=tool)) + geom_col(position = "dodge") + geom_hline(yintercept =0.90) +
    facet_wrap(~name)


df_comb <- rbind(df |> filter(tool=="deepsomatic") |> add_column(type="real"), df_simulated |> filter(tool=="deepsomatic") |> add_column(type="simulated"))
df_comb |> filter(Filter=="PASS") |> 
  ggplot(aes(x=Recall,y=Precision)) + geom_path(aes(linetype=type),linewidth=1) +
  geom_point(size=4,aes(color=purity)) + scale_color_manual(values=mypal)

df_comb |> filter(Filter=="PASS") |> pivot_longer(c(Recall,Precision)) |> 
  ggplot(aes(x=purity,y=value,fill=type)) + geom_col(position = "dodge") + geom_hline(yintercept =0.90) +
  facet_wrap(~name)


df <- df |> add_column(bc="Hom6")
df2 <- read_tsv("hcc1937_bc_hom100.tsv") |> add_column(bc="Hom100")

rbind(df,df2) |> filter(Filter=="PASS") |> pivot_longer(c(Recall,Precision)) |> 
  ggplot(aes(x=purity,y=value, fill=bc)) + geom_col(position = "dodge") + facet_grid(name~tool) + 
  geom_text(aes(label = value*100, group = bc), vjust = -0.1,position = position_dodge(width = .8), size=3) +
  labs(x="",y="",fill="Homopolymer filter")
