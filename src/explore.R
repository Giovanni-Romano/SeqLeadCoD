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
names(palette) =  sort(table(readRDS("data/rds/dataGHE.RDS")$CauseS %>% as.character()), decreasing = T) %>% names
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






# Age-year specific heterogeneity ----
m = data %>% select(c("Age", "CauseS")) %>%  
  distinct() %>% 
  group_by(Age) %>% summarise(m = n()) %>% pull(m)
Gini.max = (m-1)/m
Entropy.max = log2(m)

plt_heter = data %>% 
  group_by(Year, Age, CauseS) %>% 
  summarize(f = n()/183,
            tmp = -f*log2(f)) %>% 
  group_by(Year, Age) %>% 
  summarize(m.year = n(),
            Gini = 1 - sum(f^2),
            Entropy = sum(tmp)) %>% 
  left_join(data.frame(Age = factor(levels(data$Age), levels = levels(data$Age)),
                       Gini_max = Gini.max, Entropy_max = Entropy.max),
            by = "Age") %>% 
  mutate(Gini_normalized = Gini/Gini_max,
         Entropy_normalized = Entropy/Entropy_max,
         Nmbr_obs_causes = m.year) %>% 
  select(Year, Age, Gini_normalized, Entropy_normalized, Nmbr_obs_causes) %>% 
  ungroup() %>% 
  pivot_longer(Gini_normalized:Nmbr_obs_causes, names_to = "Index", values_to = "Heterogeneity") %>% 
  ggplot(aes(x = Age, y = Heterogeneity, col = factor(Year), group = factor(Year))) + 
  geom_line() + 
  geom_point(size = 2, alpha = 0.9) +
  geom_line(data = data.frame(Year = "Total", Age = factor(levels(data$Age), levels = levels(data$Age)), 
                              Index = "Nmbr_obs_causes", Heterogeneity = m)) +
  geom_point(data = data.frame(Year = "Total", Age = factor(levels(data$Age), levels = levels(data$Age)), 
                               Index = "Nmbr_obs_causes", Heterogeneity = m)) +
  scale_color_manual("Year", values = c(viridis::viridis(n = 22, end = 0.9), "red")) +
  guides(color = guide_legend(override.aes = list(size = 4), nrow = 2)) +
  facet_wrap(~Index, scales = "free") + 
  theme_bw() + theme(legend.position = "bottom",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(plot = plt_heter,
       filename = paste0("Age_Year_Heterogeneity_", sex,".pdf"),
       path = paste0("img/explore_", sex),
       height = 5, width = 12)


# Observed categories for each year
tmp = data %>% select(c("Age", "Year", "CauseS")) %>%  
  distinct() %>% 
  filter(Age == "5-9") %>% select(Year, CauseS) %>% 
  table
tmp.cs = colSums(tmp)
tmp[ , order(-tmp.cs)] %>% 
  pheatmap(cluster_rows = F, cluster_cols = F, 
           breaks = c(-0.5, 0.5, 1.5), 
           color = c("grey90", "black"), 
           legend_breaks = c(0, 1), legend_labels = c("Not obs.", "Observed"),
           main = "Observed causes in age class 5-9",
           height = 5, width = 12,
           filename = paste0("img/explore_", sex, "/Age_Year_ObsCauses_", sex,".pdf"))


# # NMI clustering all years ----
# hclist = list()
# entropies = matrix(NA, nrow = 50, ncol = 22)
# 
# for (y in 1:22){
#   year = 1999+y
#   df_sy = data %>% 
#     filter(Year == year) %>%  
#     unite("ID", c(CountryN, Sex)) %>% 
#     select(-c(Year, Country, CauseC, CauseT, Region, Pop2021)) %>% 
#     arrange(Age) %>% 
#     pivot_wider(names_prefix = "Age ", 
#                 names_from = Age, values_from = CauseS) %>% 
#     column_to_rownames("ID")
#   
#   dist_mat = daisy(df_sy, metric = "gower")
#   hctmp = hclust(dist_mat, method = "ward.D")
#   hclist[[y]] = hctmp
#   
#   for (k in 1:50){
#     grp = cutree(hctmp, k)
#     df_sy$cluster = grp
#     
#     # sum(sapply(tapply(df_sy, grp, daisy, metric = "gower"), mean))
#     e_tmp = df_sy %>%
#       group_by(cluster) %>%
#       summarise(across(starts_with("Age "), ~ DescTools::Entropy(table(.x)))) %>%
#       select(-cluster)
#     
#     entropies[k, y] = sum(e_tmp)/nrow(e_tmp)
#   }
# }
# 
# pdf("img/elbow_hclust.pdf", width = 12, height = 9)
# par(mfrow = c(5, 5), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))
# for (j in 1:22){
#   plot(1:50, entropies[, j], type = "l", 
#        main = paste("Year", 1999+j), 
#        xlab = "Number of clusters", ylab = "Entropy")
#   abline(v=10, col = 2, lty = 2)
# }
# dev.off()
# 
# 
# point_estimate = sapply(hclist, function(x) {
#   cutree(x, k = 10)
# })
# 
# NMI_mat = matrix(NA, nrow = 22, ncol = 22)
# for (i in 1:22){
#   for (j in 1:22){
#     NMI_mat[i, j] = aricode::NMI(point_estimate[, i], point_estimate[, j])
#   }
# }
# 
# rownames(NMI_mat) = colnames(NMI_mat) = 1999 + 1:22
# pheatmap(NMI_mat, cluster_rows = F, cluster_cols = F, breaks = seq(0, 1, by = 0.05), color = viridis::viridis(21),
#          main = "NMI between estimate from hclust with 10 clusters") %>% as.ggplot()
# ggsave("img/NMI_hclust.pdf", width = 12, height = 9)