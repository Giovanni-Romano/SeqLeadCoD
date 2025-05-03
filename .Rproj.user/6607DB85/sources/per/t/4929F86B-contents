suppressPackageStartupMessages(library(tidyverse))
library(cluster)
library(pheatmap)
library(ggplotify)
library(patchwork)

causes = readr::read_csv2("data/raw/causes.csv", 
                          col_types = cols(.default = "character"))
countries_raw = readr::read_csv2("data/raw/countries.csv", 
                                 col_types = cols(.default = "character"))
population = readr::read_csv2("data/raw/population.csv", 
                              col_types = cols(.default = "character"))
countries = left_join(countries_raw, 
                      population %>% mutate(Pop2021 = as.numeric(`2021`)) %>% select(`Country Code`, `Pop2021`),
                      by = c("Code" = "Country Code")) %>% 
  mutate(ParentCode = gsub("R", "", ParentCode) )

data = readr::read_csv2("data/raw/WHO_GHE_Top1.csv", 
                       col_types = cols(DIM_GHECAUSE_CODE = col_character())) %>% 
  select(-VAL_DTHS_RATE100K_NUMERIC) %>% 
  rename_with(~ gsub("DIM_", "", .), starts_with("DIM_")) %>% 
  rename(Country = COUNTRY_CODE, Year = YEAR_CODE,
         CauseC = GHECAUSE_CODE, CauseT = GHECAUSE_TITLE,
         Age = AGEGROUP_CODE, Sex = SEX_CODE) %>% 
  filter(!(Age %in% c("D0T27", "M1T11", "TOTAL"))) %>% 
  left_join(causes, by = c("CauseC", "CauseT")) %>%
  left_join(countries %>% 
              rename(CountryN = Title,
                     Region = ParentCode) %>% 
              select(Code, CountryN, Region, Pop2021), by = c("Country" = "Code")) %>% 
  mutate(Sex = if_else(Sex == "FEMALE", "F", "M"),
         Age = case_when(Age == "YGE_85" ~ "85+",
                         .default = gsub("Y(\\d{1,2})T(\\d{1,2})", "\\1-\\2", Age)),
         mutate(across(where(is.character), as.factor))) %>% 
  select(Year, Country, CountryN, Age, Sex, CauseC, CauseS, CauseT, Region, Pop2021)

data$Age = factor(data$Age, levels = c("0-1", "1-4", "5-9", 
                                       "10-14", "15-19", "20-24", 
                                       "25-29", "30-34", "35-39", 
                                       "40-44", "45-49", "50-54", 
                                       "55-59", "60-64", "65-69", 
                                       "70-74", "75-79", "80-84", 
                                       "85+"))



# Clustering single year ----
year = 2000
df_sy = data %>% 
  filter(Year == year) %>%  
  unite("ID", c(CountryN, Sex)) %>% 
  select(-c(Year, Country, CauseC, CauseT, Region, Pop2021)) %>% 
  arrange(Age) %>% 
  pivot_wider(names_prefix = "Age ", 
              names_from = Age, values_from = CauseS) %>% 
  column_to_rownames("ID")

dist_mat = daisy(df_sy, metric = "gower")
hc = hclust(dist_mat, method = "ward.D")

ncauses = length(unique(data %>% select(CauseS) %>% unlist))
palette = c(pals::alphabet(26)[-(4:5)], pals::alphabet2(26), pals::glasbey(32))[1:ncauses]
names(palette) =  sort(table(data$CauseS %>% as.character()), decreasing = T) %>% names


for (yyyy in 2000:2021){
plt = data %>% 
  filter(Year == yyyy) %>%
  unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
  mutate(ID = factor(ID, levels = rownames(df_sy)[hc$order]), .after = ID) %>% 
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
                    mutate(ID = factor(ID, levels = rownames(df_sy)[hc$order]), .after = ID) %>% 
                    select(ID, Sex, Region)) %>% pivot_longer(cols = -ID, names_to = "x", values_to = "fill") %>% 
  mutate(x = factor(x, levels = c("Sex", "Region")))

colors_aux = c("lightblue", "pink", "#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")
names(colors_aux) = c("M", "F", "AF", "AM", "EM", "EU", "SEA", "WP")

plt_aux = ggplot(data = data_aux) + 
  geom_tile(aes(x = x, y = ID, fill = fill)) +
  scale_fill_manual("Region", values = colors_aux) +
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

ggsave(plot = plt_out,
       filename = paste0("heatmap_", yyyy, ".pdf"),
       path = "img/heatmaps",
       height = 40, width = 16)
}


# Clustering single age_classes ----
ageclass = levels(data$Age)[19]
df_sa = data %>% 
  filter(Age == ageclass) %>%  
  unite("ID", c(CountryN, Sex)) %>% 
  select(-c(Age, Country, CauseC, CauseT, Region, Pop2021)) %>% 
  arrange(Year) %>% 
  pivot_wider(names_prefix = "Year ", 
              names_from = Year, values_from = CauseS) %>% 
  column_to_rownames("ID")

dist_mat.ages = daisy(df_sa, metric = "gower")
hc.ages = hclust(dist_mat.ages, method = "ward.D")

for (aaaa in levels(data$Age)){
  plt = data %>% 
    filter(Age == aaaa) %>%
    unite("ID", c(CountryN, Sex), remove = FALSE) %>% 
    mutate(ID = factor(ID, levels = rownames(df_sa)[hc.ages$order]), .after = ID) %>% 
    ggplot(aes(x = Year, y = ID, fill = CauseS)) +
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
                      mutate(ID = factor(ID, levels = rownames(df_sa)[hc.ages$order]), .after = ID) %>% 
                      select(ID, Sex, Region)) %>% pivot_longer(cols = -ID, names_to = "x", values_to = "fill") %>% 
    mutate(x = factor(x, levels = c("Sex", "Region")))
  
  plt_aux = ggplot(data = data_aux) + 
    geom_tile(aes(x = x, y = ID, fill = fill)) +
    scale_fill_manual("Region", values = colors_aux) +
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
  
  ggsave(plot = plt_out,
         filename = paste0("heatmap_", aaaa, ".pdf"),
         path = "img/heatmaps_ages_orderfixed",
         height = 40, width = 24)
}


# Plots grids USA male and female ----
plt_USA.M = data %>% 
  filter(Country == "USA", Sex == "M") %>% 
  ggplot(aes(x = Age, y = Year, fill = CauseS)) +
  geom_raster() + 
  geom_text(aes(label =CauseS), size = 2) +
  scale_fill_manual(values = palette) +
  ggtitle("Male population in USA") + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
plt_USA.F = data %>% 
  filter(Country == "USA", Sex == "F") %>% 
  ggplot(aes(x = Age, y = Year, fill = CauseS)) +
  geom_raster() + 
  geom_text(aes(label =CauseS), size = 2) +
  scale_fill_manual(values = palette) +
  ggtitle("Male population in USA") + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = plt_USA.M,
       filename = "heatmap_USAMale.pdf",
       path = "img/",
       height = 6, width = 14)
ggsave(plot = plt_USA.F,
       filename = "heatmap_USAFemale.pdf",
       path = "img/",
       height = 6, width = 14)

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