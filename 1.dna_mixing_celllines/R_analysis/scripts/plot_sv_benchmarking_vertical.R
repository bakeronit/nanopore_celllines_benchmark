library(here)

## after loading df_sv in 3.benchmark_sv_calling.Rmd

p1v <- df_sv |> filter(cellline=="COLO829")|>  mutate(ID=factor(ID,levels=id_levels), purity=factor(purity,levels=seq(10,100,10) |> as.character())) |> 
  ggplot(aes(x=purity,y=ID,fill=svtype, alpha=supp)) + 
  geom_tile(color="black") + 
  facet_grid2(svtype~tool, scales = "free", space = "free",
              strip = strip_themed(background_y = elem_list_rect(fill=c("#BC4749","#F4ACB7","#588157" ,"#FF9F1C")),
                                   background_x = elem_list_rect(fill="white"))
  ) +
  scale_alpha_manual(values = c(0.1, 1), labels = c("Missed", "Detected")) + 
  scale_fill_manual(values = sv_palette) + theme_bw(base_size = 12) +
  theme(legend.position="bottom",
        legend.direction = "horizontal",     
        legend.margin = margin(0, 0, 0, 0), axis.text.y.left  = element_blank(), axis.ticks.y = element_blank())+  # Remove margins around the legend
  #scale_x_discrete(labels = function(x) str_remove(x, "0_0_truthset_")) +
  labs(y="SV gold standard", x="Purity", alpha="", fill="SV type")

 p1v
 p1v <- p1v + theme(legend.position = "none")
 ggsave(here("figures/sv_heatmap_colo829_vertical.png"), p1v,height = 7.6,width = 8.4) 
 
 
 p2v <- df_sv |> filter(cellline=="HCC1937")|>  mutate(ID=factor(ID,levels=id_levels), purity=factor(purity,levels=seq(10,100,10) |> as.character())) |> 
   ggplot(aes(x=purity,y=ID,fill=svtype, alpha=supp)) + 
   geom_tile(color="black") + 
   facet_grid2(svtype~tool, scales = "free", space = "free",
               strip = strip_themed(background_y = elem_list_rect(fill=c("#BC4749","#F4ACB7","#588157" ,"#FF9F1C")),
                                    background_x = elem_list_rect(fill="white"))
   ) +
   scale_alpha_manual(values = c(0.1, 1), labels = c("Missed", "Detected")) + 
   scale_fill_manual(values = sv_palette) + theme_bw(base_size = 12) +
   theme(legend.position="bottom",
         legend.direction = "horizontal",     
         legend.margin = margin(0, 0, 0, 0), axis.text.y.left  = element_blank(), axis.ticks.y = element_blank())+  # Remove margins around the legend
   #scale_x_discrete(labels = function(x) str_remove(x, "0_0_truthset_")) +
   labs(y="SV gold standard", x="Purity", alpha="", fill="SV type")
 
 p2v
 p2v <- p2v + theme(legend.position = "none")
 ggsave(here("figures/sv_heatmap_hcc1937_vertical.png"), p2v, height = 12.6,width = 8.4)
 