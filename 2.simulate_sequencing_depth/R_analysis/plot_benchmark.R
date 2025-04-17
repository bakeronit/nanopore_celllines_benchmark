library(tidyverse)
bc_path <- "/working/lab_nicw/jiaZ/bioprojects/nanopore_celllines_benchmark/2.simulate_sequencing_depth/analysis/benchmark/snvs/somatic/clairS"
df <- list.files(bc_path, pattern = "summary.txt",recursive = TRUE, full.names = TRUE) |> 
  set_names(\(x) str_split_i(x,"/",13)) |> map(\(file) read.table(file,header=TRUE)) |> list_rbind(names_to = "prefix") |> 
  separate_wider_delim(prefix,delim = ".",names = c("sample_t","depth_t","sample_n","depth_n"))


df |> filter(Filter=="PASS") |> unite(col = "sequencing_plan", c(depth_t,depth_n)) |> select(sample_t,sequencing_plan,Precision,Recall) |> 
  pivot_longer(c(Precision,Recall)) |> 
  mutate(sample_t=factor(sample_t, levels=c("COLO829_40","COLO829_60","COLO829_80","COLO829"),labels=c("40%","60%","80%","100%"))) |> 
  ggplot(aes(x=sample_t,y=value, fill=name)) + 
  geom_col(position="dodge") + 
  geom_text(aes(label = value*100, group = name), vjust = -0.1,position = position_dodge(width = .8), size=3) + 
  facet_wrap(~sequencing_plan, nrow = 3) + labs(x="",y="",fill="")
