suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(ggalluvial)
  library(patchwork)
})
source("src/tHMM_utils.R")


# Simulate data ----
out = simdata(seed = 1)
C = out$C
Theta = out$Theta
Y = out$Y
nT = out$info$nT

# Palette ----
palette = palette = c("#9DCC00", "#808080",
                      "#005C31",
                      "#993F00", "#0075DC",
                      "#F0A0FF", "#5EF1F2",
                      "#FFCC99", "#782AB6", 
                      "#FFA405")
names(palette) = 1:10


# Plots ----
Cflow = C %>% melt(value.name = "Cluster") %>%
  mutate(freq = 1) %>% 
  mutate(Cluster = factor(Cluster, levels = 1:max(C, na.rm = T)))

plt_flow = Cflow %>% ggplot(aes(x = Year, stratum = Cluster, alluvium = Subject,
                                y = freq,
                                fill = Cluster)) +
  geom_flow() +
  geom_stratum(color = "white") + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3) + 
  coord_cartesian(expand = TRUE, xlim = c(0.95, 8.05)) +
  scale_x_continuous(breaks = 1:nT) +
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

plt_th = Theta %>% melt(value.name = "Cause") %>%
  mutate(Cause = factor(Cause, levels = 1:max(Theta, na.rm = T))) %>%
  ggplot(aes(x = Age, y = -Cluster, fill = Cause)) +
  geom_tile() +
  geom_text(aes(label = Cause), size = 3, col = "white") +
  scale_fill_manual(values = palette, breaks = names(palette), na.value = "transparent") +
  facet_wrap(~ Year, scales = "free_y", nrow = 1) +
  labs(x = "Age", y = "Clusters", fill = "Cause") + 
  theme(axis.text.y =  element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(margin = margin(2, 2, 2, 2)))

plt_param = (plt_flow/plt_th) + plot_layout (heights = c(2, 1))


list_plt = lapply(1:(nT), function(t) 
  Y %>% 
    reshape2::melt(value.name = "value") %>% 
    mutate(value = factor(value)) %>% 
    left_join(reshape2::melt(C) %>% rename(Cluster = value), by = c("Subject", "Year")) %>%
    unite("ID", c("Subject", "Year"), remove = F) %>%
    filter(Year == t) %>% 
    ggplot(aes(x = Age, y = ID, fill = value)) +
    geom_tile() +
    scale_fill_manual(values = palette) +
    # ggh4x::facet_grid2(rows = vars(Cluster), cols = vars(Year), scales = "free", independent = "all")
    facet_grid(rows=vars(Cluster), cols = vars(Year), scales = "free_y", space = "free")+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(margin = margin(2, 2, 2, 2)),
          panel.spacing = unit(0.1, "lines"),  # smaller space between plots
          plot.margin = margin(2, 2, 2, 2)     # tighten outer margin per plot
    ))

library(gridExtra)
library(grid)
# Arrange the 8 plots (2 rows × 4 columns)
plot_matrix = arrangeGrob(
  grobs = list_plt,
  nrow = 1,
  ncol = 8
)

# Create common axis titles as grobs
top_label = textGrob("Age", gp = gpar(fontsize = 14))
bottom_label = textGrob("Year", gp = gpar(fontsize = 14))
right_label = textGrob("Cluster", rot = 90, gp = gpar(fontsize = 14))
left_label = textGrob("Population", rot = -90, gp = gpar(fontsize = 14))

# Combine all into a master layout
plt_data = grid.arrange(
  arrangeGrob(
    left = left_label,
    right = right_label,
    top = top_label,
    bottom = bottom_label,
    plot_matrix
  )
)


IMGFOLDER = "img/simstudy4"
if (!dir.exists(IMGFOLDER)) {
  dir.create(IMGFOLDER, recursive = TRUE)
}

ggsave(path = IMGFOLDER, 
       filename = "tHMM_simplot_param.pdf",
       plot = plt_param,
       width = 10, height = 6)

ggsave(path = IMGFOLDER,
       filename = "tHMM_simplot_data.pdf",
       plot = plt_data,
       width = 10, height = 6)
  