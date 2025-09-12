# Import pkgs ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(patchwork)
})

# Import souce code
source("src/utils.R")

# Load data ----
sex = "female"
if (sex == "female"){
  data = readRDS("data/rds/GHEdf_female.RDS")
} else if (sex == "male"){
  data = readRDS("data/rds/GHEdf_male.RDS")
} else {
  stop("Please select a value for sex among (female, male)")
}

causes = readr::read_csv2("data/raw/causes.csv") %>% select(CauseS) %>% unique() %>% unlist() %>% unname() %>% sort() 
ncauses = length(causes)
palette = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(32))[1:ncauses]
names(palette) =  sort(table(data$CauseS %>% as.character()), decreasing = T) %>% names
nregions = levels(data$Region) %>% length
palette_regions = Polychrome::glasbey.colors(nregions + 1)[-1]
names(palette_regions) = levels(data$Region)
colortext_regions = c("white", "black", "black", "white", "black", "white", "black", "black")
names(colortext_regions) = levels(data$Region)


for (yyyy in 2000:2021){
  
  dist_mat = data %>% 
    filter(Year == yyyy) %>%  
    select(-c(Year, Sex, Country, CauseC, CauseT, Region)) %>% 
    arrange(Age) %>% 
    pivot_wider(names_prefix = "Age ", 
                names_from = Age, values_from = CauseS) %>% 
    column_to_rownames("CountryN") %>% 
    hamming_matrix()
  
  hc = hclust(dist_mat, method = "ward.D")
  
  plt = data %>% 
    filter(Year == yyyy) %>%
    mutate(CountryN = factor(CountryN, levels = rownames(df_sy)[hc$order])) %>% 
    ggplot(aes(y = Age, x = CountryN, fill = CauseS)) +
    geom_raster() + 
    geom_text(aes(label = CauseS), size = 2, angle = 90) +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())
  
  data_aux = unique(data %>% select(CountryN, Region)) %>% 
    mutate(CountryN = factor(CountryN, levels = rownames(df_sy)[hc$order])) %>% 
    cbind(., x = "Region")
  
  
  plt_aux = ggplot(data = data_aux,
                   mapping = aes(x = x, y = CountryN, fill = Region)) + 
    geom_tile() +
    scale_fill_manual("Region", values = palette_regions) +
    scale_color_manual(values = colortext_regions) +
    geom_text(aes(label = Region, color = Region), size = 2, angle = 90) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white"),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title.x = element_blank()) +
    coord_flip()
  
  plt_out = plt / plt_aux + plot_layout(heights = c(20, 1))
  
  ggsave(plot = plt_out,
         filename = paste0("heatmap_", yyyy, ".pdf"),
         path = paste0("img/explore_", sex, "/heatmaps"),
         height = 12, width = 20)
}










# NMI clustering all years ----
hclist = list()
entropies = matrix(NA, nrow = 50, ncol = 22)

for (y in 1:22){
  year = 1999+y
  df_sy = data %>% 
    filter(Year == year) %>%  
    unite("ID", c(CountryN, Sex)) %>% 
    select(-c(Year, Country, CauseC, CauseT, Region, Pop2021)) %>% 
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
    
    # sum(sapply(tapply(df_sy, grp, daisy, metric = "gower"), mean))
    e_tmp = df_sy %>%
      group_by(cluster) %>%
      summarise(across(starts_with("Age "), ~ DescTools::Entropy(table(.x)))) %>%
      select(-cluster)
    
    entropies[k, y] = sum(e_tmp)/nrow(e_tmp)
  }
}

pdf("img/elbow_hclust.pdf", width = 12, height = 9)
par(mfrow = c(5, 5), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))
for (j in 1:22){
  plot(1:50, entropies[, j], type = "l", 
       main = paste("Year", 1999+j), 
       xlab = "Number of clusters", ylab = "Entropy")
  abline(v=10, col = 2, lty = 2)
}
dev.off()


point_estimate = sapply(hclist, function(x) {
  cutree(x, k = 10)
})

NMI_mat = matrix(NA, nrow = 22, ncol = 22)
for (i in 1:22){
  for (j in 1:22){
    NMI_mat[i, j] = aricode::NMI(point_estimate[, i], point_estimate[, j])
  }
}

rownames(NMI_mat) = colnames(NMI_mat) = 1999 + 1:22
pheatmap(NMI_mat, cluster_rows = F, cluster_cols = F, breaks = seq(0, 1, by = 0.05), color = viridis::viridis(21),
         main = "NMI between estimate from hclust with 10 clusters") %>% as.ggplot()
ggsave("img/NMI_hclust.pdf", width = 12, height = 9)