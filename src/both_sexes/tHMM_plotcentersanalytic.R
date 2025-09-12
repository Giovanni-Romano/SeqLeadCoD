ncl_df = ghe_PPE %>% 
  group_by(Year, Cl_new) %>% 
  summarise(n_cl = n()/19)

ncauses_df = ghe_PPE %>% 
  group_by(Year, Cl_new, Age, CauseS) %>% 
  summarise(n_cause = n())

center_df.analytic = left_join(ncauses_df, ncl_df, by = c("Year", "Cl_new")) %>% 
  mutate(p_cause = n_cause / n_cl) %>%
  group_by(Year, Cl_new, Age) %>% 
  slice_max(order_by = p_cause, n = 1) %>% 
  ungroup


center_df.both = full_join(center_df.analytic %>% 
                             rename(CauseAnalytic = CauseS), 
                           center_df %>% rename(CauseSample = Cause), 
                           by = c("Year", "Cl_new", "Age")) 

mean(center_df.both$CauseAnalytic == center_df.both$CauseSample)

all_clusters = unique(c(PPE))

center_df.complete = center_df.both %>% 
  mutate(order = as.integer(as.character(order))) %>% 
  complete(Year, Age, Cl = all_clusters) %>%
  # Impute order also for "fake" clusters
  group_by(Year, Age) %>% 
  mutate(order = ifelse(is.na(order),
                        max(order, na.rm = TRUE) + dense_rank(Cl[is.na(order)]),
                        order)) %>%
  ungroup() %>% 
  mutate(order = factor(order, levels = 484:1))

plt_centers.analytic = center_df.complete %>% 
  ggplot() +
  geom_raster(aes(x = as.integer(Age), y = order, fill = CauseAnalytic)) +
  scale_fill_manual(values = palette, na.value = "transparent",
                    breaks = sort(names(palette))) +
  guides(fill = guide_legend(nrow = 6, override.aes = list(color = "black"))) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), 
                     labels = ages[seq(1, nages, by = 3)]) +
  facet_wrap(~ Year, ncol = 5, scales = "free_y") +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),  
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    text = element_text(size = 15))

ggsave(plt_centers.analytic, 
       filename = "img/tHMM/Gnedin20019/centers_analytic.pdf",
       width = 20*0.7, height = 20*0.7)
