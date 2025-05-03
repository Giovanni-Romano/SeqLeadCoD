library(tidyverse)
library(patchwork)
library(magick)

# Create a list to store plots
plot_list <- list()

# Loop over years
for (yyyy in 2000:2021) {
  df <- data %>% filter(Year == yyyy)
  
  plt <- data %>%
    filter(Year == yyyy) %>%
    unite("ID", c(CountryN, Sex), remove = FALSE) %>%
    mutate(ID = factor(ID, levels = rownames(df_sy)[hc$order]), .after = ID) %>%
    ggplot(aes(x = Age, y = ID, fill = CauseS)) +
    geom_raster() +
    scale_y_discrete(name = "Country+Sex") +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.ticks = element_blank(),
          axis.text = element_blank()) +
    ggtitle(as.character(yyyy))  # Optional: show the year as a frame title
  
  plot_list[[as.character(yyyy)]] <- plt
}

# Convert plots to images
img_frames <- lapply(plot_list, function(p) {
  img <- magick::image_graph(width = 1600, height = 1600, res = 100)
  print(p)
  dev.off()
  img
})

# Create and write animated gif
animation <- image_animate(image_join(img_frames), fps = 1)
image_write(animation, path = "img/heatmaps/heatmap_animation.gif")
