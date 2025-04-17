a <-fp_mut |> rownames_to_column("sample_id") |> pivot_longer(-sample_id) |> 
  separate_wider_delim(sample_id, names= c("sample","purity","tool"), delim = "@") |> 
  mutate(purity=factor(purity, levels=seq(100,10,-10))) |> 
  mutate(c1=str_split_i(name, "\\[", i=1), c2=str_split_i(name,"\\]",i=2), mut=sub(".*\\[(.*?)\\].*", "\\1", name)) |> 
  unite("context",c(c1,c2), sep = "") |> 
  group_by(sample, purity, tool) |> mutate(sum=sum(value), prop=value/sum)


gc <- gc_colo829 |> rownames_to_column("sample") |> pivot_longer(-sample) |> 
  mutate(gc=value/sum(value)) |> select(-c(sample,value))

a |> left_join(gc, by="name") |> filter(startsWith(sample,"COLO829")) |> 
  ggplot(aes(x=log(gc), y=log(prop), color=mut)) +
  geom_point(alpha=.6,size=2) + ggrepel::geom_text_repel(aes(label=name)) +
  facet_grid(tool~purity, scales = "free") + scale_color_manual(values = COLORS6)
