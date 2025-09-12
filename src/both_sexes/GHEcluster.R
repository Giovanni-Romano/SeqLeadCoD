suppressPackageStartupMessages(library(tidyverse))
library(cluster)
library(pheatmap)
library(ggplotify)
library(patchwork)

data = readRDS("data/rds/dataGHE.RDS")

# Palette ----
ncauses = length(unique(data %>% select(CauseS) %>% unlist))
palette = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(32))[1:ncauses]
names(palette) =  sort(table(data$CauseS %>% as.character()), decreasing = T) %>% names

# NMI clustering all years ----
hclist = list()
entropies = matrix(NA, nrow = 50, ncol = 22)

for (y in 1:22){
  year = 1999+y
  df_sy = data %>% 
    filter(Year == year) %>%  
    unite("ID", c(CountryN, Sex)) %>% 
    select(-c(Year, Country, CauseC, CauseT, Region)) %>% 
    arrange(Age) %>% 
    pivot_wider(names_prefix = "Age ", 
                names_from = Age, values_from = CauseS) %>% 
    column_to_rownames("ID")
  
  dist_mat = daisy(df_sy, metric = "gower")
  hctmp = hclust(dist_mat, method = "ward.D")
  hclist[[y]] = hctmp
  
  for (k in 1:50){
    grp = cutree(hctmp, k)
    df_sy$cluster = grp
    
    e_tmp = df_sy %>%
      group_by(cluster) %>%
      summarise(across(starts_with("Age "), ~ DescTools::Entropy(table(.x)))) %>%
      select(-cluster)
    
    entropies[k, y] = sum(e_tmp)/nrow(e_tmp)
  }
}

ncl_hc = 15
par(mfrow = c(5, 5), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))
for (j in 1:22){
  plot(1:50, entropies[, j], type = "l", 
       main = paste("Year", 1999+j), 
       xlab = "Number of clusters", ylab = "Entropy")
  abline(v=ncl_hc, col = 2, lty = 2)
}


point_estimate = sapply(hclist, function(x) {
  tmp = cutree(x, k = ncl_hc)
  tmp = tmp[order(names(tmp))]
  tmp
})
colnames(point_estimate) = 2000:2021

NMI_mat = matrix(NA, nrow = 22, ncol = 22)
for (i in 1:22){
  for (j in 1:22){
    NMI_mat[i, j] = aricode::NMI(point_estimate[, i], point_estimate[, j])
  }
}

rownames(NMI_mat) = colnames(NMI_mat) = 1999 + 1:22
pheatmap(NMI_mat, cluster_rows = F, cluster_cols = F, breaks = seq(0, 1, by = 0.05), color = viridis::viridis(21),
         main = paste0("NMI between estimate from hclust with ", ncl_hc, " clusters")) %>% as.ggplot()


# Plots with column for cluster estimates (rows ordered according to partition in 2000) ----
point_estimate.df = point_estimate %>% as.data.frame() %>% rownames_to_column("ID") %>% 
  pivot_longer(-ID, names_to = "Year", values_to = "HC") %>% 
  mutate(Year = as.numeric(Year), HC = factor(as.character(HC), levels = sort(unique(HC))))

hc2000 = hclist[[1]]
countries2000 = hc2000$labels


plt_order2000 = list()
for (y in 1:22){
  yyyy = 1999+y
  plt = data %>% 
    filter(Year == yyyy) %>%
    unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
    mutate(ID = factor(ID, levels = countries2000[hc2000$order]), .after = ID) %>% 
    ggplot(aes(x = Age, y = ID, fill = CauseS)) +
    geom_raster() + 
    geom_text(aes(label = CauseS), size = 2) +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ppeyyyy = point_estimate.df %>% filter(Year == yyyy) %>% select(-Year)
  data_aux = unique(data %>% unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
                      select(ID, Sex, Region)) %>% 
    left_join(ppeyyyy, by = "ID") %>%
    mutate(ID = factor(ID, levels = countries2000[hc2000$order]), .after = ID) %>% 
    pivot_longer(cols = -ID, names_to = "x", values_to = "fill") %>% 
    mutate(x = factor(x, levels = c("Sex", "Region", "HC"), labels = c("Sex", "Region", "Cluster")))
  
  colors_aux = c(c("lightblue", "pink"), 
                 c("#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C"),
                 rainbow(length(unique(ppeyyyy$HC))))
  names(colors_aux) = c(c("M", "F"), 
                        c("AF", "AM", "EM", "EU", "SEA", "WP"),
                        as.character(1:ncl_hc))
  
  plt_aux = ggplot(data = data_aux) + 
    geom_tile(aes(x = x, y = ID, fill = fill)) +
    scale_fill_manual(values = colors_aux) +
    geom_text(aes(x = x, y = ID, label = fill), size = 2) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  plt_out = plt_aux + plt + plot_layout(widths = c(1, 20))
  
  plt_order2000[[y]] = plt_out
}


# Plots with column for cluster estimates (rows ordered according to partition of corresponding year) ----
plt_orderfree = list()
for (y in 1:22){
  yyyy = 1999+y
  
  ppeyyyy = point_estimate.df %>% filter(Year == yyyy) %>% select(-Year)
  hcyyyy = hclist[[y]]
  
  plt = data %>% 
    filter(Year == yyyy) %>%
    unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
    mutate(ID = factor(ID, levels = hcyyyy$labels[hcyyyy$order]), .after = ID) %>% 
    ggplot(aes(x = Age, y = ID, fill = CauseS)) +
    geom_raster() + 
    geom_text(aes(label = CauseS), size = 2) +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  data_aux = unique(data %>% unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
                      select(ID, Sex, Region)) %>% 
    left_join(ppeyyyy, by = "ID") %>%
    mutate(ID = factor(ID, levels = hcyyyy$labels[hcyyyy$order]), .after = ID) %>% 
    pivot_longer(cols = -ID, names_to = "x", values_to = "fill") %>% 
    mutate(x = factor(x, levels = c("Sex", "Region", "HC"), labels = c("Sex", "Region", "Cluster")))
  
  colors_aux = c(c("lightblue", "pink"), 
                 c("#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C"),
                 rainbow(length(unique(ppeyyyy$HC))))
  names(colors_aux) = c(c("M", "F"), 
                        c("AF", "AM", "EM", "EU", "SEA", "WP"),
                        as.character(1:10))
  
  plt_aux = ggplot(data = data_aux) + 
    geom_tile(aes(x = x, y = ID, fill = fill)) +
    scale_fill_manual(values = colors_aux) +
    geom_text(aes(x = x, y = ID, label = fill), size = 2) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())
  
  plt_out = plt_aux + plt + plot_layout(widths = c(1, 20))
  
  plt_orderfree[[y]] = plt_out
}
