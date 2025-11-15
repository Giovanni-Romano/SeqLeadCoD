options(warn = 1)
# Packages ----
suppressPackageStartupMessages(library(tidyverse))
library(abind)
library(ggh4x)
library(grid)
library(gridExtra)
# pkgs for plotting maps
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
source("src/utils.R")

# Settings ----
theme_set(theme_bw())
IMGFOLDER = paste0("img/tHMM/")

# Load results ----
res_male <- readRDS(paste0("output/tHMM/5yrs/PT/male/res_Gnedin20146.RDS"))
res_female <- readRDS(paste0("output/tHMM/5yrs/PT/female/res_Gnedin20019.RDS"))

PPE_male <- readRDS(paste0("output/tHMM/5yrs/PT/male/PPE_all.RDS"))[["20146"]]
PPE_female <- readRDS(paste0("output/tHMM/5yrs/PT/female/PPE_all.RDS"))[["20019"]]

# Load datasets
ghe_F = readRDS("data/rds/GHEdf_female_5yrs.RDS")
ghe_M = readRDS("data/rds/GHEdf_male_5yrs.RDS")
ghe = rbind(ghe_F, ghe_M)

## Define dimensions ----
ages = levels(ghe$Age)
nages = length(ages)
years = seq(2000, 2020, by = 5)
nyears = length(years)
country_codes = ghe %>% pull(Country) %>% unique %>% as.character() %>% sort
ncountries = length(country_codes)

# Palette ----
palette_regions = c(
  "Western Offshoots" = "#EF476F",
  "Latin America & Caribbean" = "#F78C6B",
  "Post-Sovietic" = "#6a994e",
  "Europe (Non P-S)" = "#a7c957",
  "Middle East & North Africa" = "#a68a64",
  "Sub-Saharan Africa" = "#d4a017",
  "South Asia" = "#8ecae6",
  "East Asia & Pacific" = "#118AB2"
)


# Pairwise hamming distances between rows ----
HSM_M = lapply(years, function(yyyy) {
  HDM = ghe_M %>% filter(Year == yyyy) %>% 
    select(Country, Age, CauseS) %>% 
    arrange(Country, Age) %>% 
    pivot_wider(id_cols = Country, names_from = Age, values_from = CauseS) %>%
    column_to_rownames("Country") %>% 
    hamming_matrix() %>% as.matrix()
  HSM = (19 - HDM)/19}
) %>% abind(., along = 3)

HSM_F = lapply(years, function(yyyy) {
  HDM = ghe_F %>% filter(Year == yyyy) %>% 
    select(Country, Age, CauseS) %>% 
    arrange(Country, Age) %>% 
    pivot_wider(id_cols = Country, names_from = Age, values_from = CauseS) %>% 
    column_to_rownames("Country") %>% 
    hamming_matrix() %>% as.matrix()
  HSM = (19 - HDM)/19}
) %>% abind(., along = 3)

dimnames(HSM_M) = dimnames(HSM_F) = list("Country1" = country_codes,
                                         "Country2" = country_codes,
                                         "Year" = years)


# HSM plot ----
## Males ---- 
ph2000_M = pheatmap::pheatmap(HSM_M[ , , "2000"], 
                              clustering_distance_rows = as.dist(1-HSM_M[ , , "2000"]),
                              clustering_distance_cols = as.dist(1-HSM_M[ , , "2000"]),
                              clustering_method = "ward.D")

HSM_M.dfplot = HSM_M %>% reshape2::melt() %>% 
  left_join(ghe_M %>% select(Year, Country, Region) %>% unique() %>% rename(Region1 = Region), by = c("Year", "Country1" = "Country")) %>%
  left_join(ghe_M %>% select(Year, Country, Region) %>% unique() %>% rename(Region2 = Region), by = c("Year", "Country2" = "Country")) %>%
  mutate(Country1 = factor(Country1, levels = ph2000_M$tree_row$labels[ph2000_M$tree_row$order], ordered = T),
         Country2 = factor(Country2, levels =  ph2000_M$tree_row$labels[ph2000_M$tree_row$order], ordered = T))

listplot_HSM_M = lapply(years, function(y)
  HSM_M.dfplot %>%
    filter(Year == y) %>% 
    ggplot() +
    geom_raster(aes(x = fct_rev(Country1), y = Country2, fill = value)) +
    scale_fill_gradient(low = "#edf6f9", high = "#264653", limits = c(0, 1), na.value = "red") +
    facet_nested(rows = vars(Region2), cols = vars(Year, Region1), scales = "free", space = "free",
                 strip = strip_nested(text_x = append(list(element_text(size = 16)), replicate(8, element_text(size = 8), simplify = F)),
                                      text_y = replicate(8, element_text(size = 8), simplify = F))
    ) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          strip.text = element_text(size = 8),
          legend.position = "none",
          panel.spacing = unit(0.1, "lines"))
)

width = 18; height = 4; zoom = 1
ggsave("Intro_HSM_M.pdf",
       path = IMGFOLDER,
       arrangeGrob(grobs = listplot_HSM_M,
                   top = grid::textGrob(
                     "Hamming similarity matrices for male populations",
                     x = 0.0055,            # left alignment in the figure
                     hjust = 0,        # justification to the left
                     gp = grid::gpar(fontsize = 20) # optional styling
                   ),
                   nrow = 1),
       width = width * zoom, height = height * zoom)


## Females ----
ph2000_F = pheatmap::pheatmap(HSM_F[ , , "2000"], 
                              clustering_distance_rows = as.dist(1-HSM_F[ , , "2000"]),
                              clustering_distance_cols = as.dist(1-HSM_F[ , , "2000"]),
                              clustering_method = "ward.D")

HSM_F.dfplot = HSM_F %>% reshape2::melt() %>% 
  left_join(ghe_F %>% select(Year, Country, Region) %>% unique() %>% rename(Region1 = Region), by = c("Year", "Country1" = "Country")) %>%
  left_join(ghe_F %>% select(Year, Country, Region) %>% unique() %>% rename(Region2 = Region), by = c("Year", "Country2" = "Country")) %>%
  mutate(Country1 = factor(Country1, levels = ph2000_F$tree_row$labels[ph2000_F$tree_row$order], ordered = T),
         Country2 = factor(Country2, levels =  ph2000_F$tree_row$labels[ph2000_F$tree_row$order], ordered = T))

listplot_HSM_F = lapply(years, function(y)
  HSM_F.dfplot %>%
    filter(Year == y) %>% 
    ggplot() +
    geom_raster(aes(x = fct_rev(Country1), y = Country2, fill = value)) +
    scale_fill_gradient(low = "#edf6f9", high = "#264653", limits = c(0, 1), na.value = "red") +
    facet_nested(rows = vars(Region2), cols = vars(Year, Region1), scales = "free", space = "free",
                 strip = strip_nested(text_x = append(list(element_text(size = 16)), replicate(8, element_text(size = 8), simplify = F)),
                                      text_y = replicate(8, element_text(size = 8), simplify = F))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          panel.spacing = unit(0.1, "lines"))
)

width = 18; height = 4; zoom = 1
ggsave("Intro_HSM_F.pdf",
       path = IMGFOLDER,
       arrangeGrob(grobs = listplot_HSM_F,
                   top = grid::textGrob(
                     "Hamming similarity matrices for female populations",
                     x = 0.0055,            # left alignment in the figure
                     hjust = 0,        # justification to the left
                     gp = grid::gpar(fontsize = 20) # optional styling
                   ),
                   nrow = 1),
       width = width * zoom, height = height * zoom)


# Plot map ----
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  left_join(ghe_F %>% select(Country, RegionN) %>% distinct(),
            by = c("iso_a3_eh" = "Country"))

palette_regions = c(
  "Western Offshoots" = "#EF476F",
  "Latin America & Caribbean" = "#F78C6B",
  "Post-Sovietic" = "#6a994e",
  "Europe (Non P-S)" = "#a7c957",
  "Middle East & North Africa" = "#a68a64",
  "Sub-Saharan Africa" = "#d4a017",
  "South Asia" = "#8ecae6",
  "East Asia & Pacific" = "#118AB2"
)

plt_map =
  ggplot(world) +
  geom_sf(aes(fill = factor(RegionN))) +
  scale_fill_manual("Regions",
                    values = palette_regions, na.value = "white",
                    labels = function(x) str_wrap(x, width = 15),
                    breaks = levels(ghe_F$RegionN)) +
  coord_sf(crs = "+proj=robin") +
  theme_bw() +
  theme(legend.position = "right",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "grey95"),
        panel.grid = element_line(colour = "white"),
        legend.key.spacing.y = unit(3, "pt"),
        legend.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.1),
        legend.text = element_text(size = 10))

{
  width = 7.5; height = 3; zoom = 1
  filename = "Intro_map.pdf"
  ggsave(filename = filename,
         plot = plt_map,
         path = IMGFOLDER,
         width = width*zoom, height = height*zoom)
  }


# PPE similarities across years ----
## Alpha ----
alpha_male <- res_male$traces[[1]]$alpha[-1, -(1:5000)]
alpha_female <- res_female$traces[[1]]$alpha[-1, -(1:5000)]

alpha_male.est <- rbind(
  Mean = rowMeans(alpha_male),
  apply(alpha_male, 1, HDInterval::hdi)
) %>% t()

alpha_female.est <- rbind(
  Mean = rowMeans(alpha_female),
  apply(alpha_female, 1, HDInterval::hdi)
) %>% t()

dimnames(alpha_male.est) <- dimnames(alpha_female.est) <- list(
  "Transition" = paste(seq(2000, 2015, by = 5),
                       seq(2005, 2020, by = 5),
                       sep = "-"
  ),
  "Quantity" = c("Mean", "HDI_L", "HDI_U")
)

alpha_plot <- rbind(
  alpha_male.est %>% as.data.frame() %>% rownames_to_column("Transition") %>%
    mutate(Sex = "M"),
  alpha_female.est %>% as.data.frame() %>% rownames_to_column("Transition") %>%
    mutate(Sex = "F")
) %>%
  mutate(Quantity = "Smoothness parameter alpha") %>%
  ggplot(aes(x = Transition, colour = Sex)) +
  geom_point(aes(y = Mean), size = 4, position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = HDI_L, ymax = HDI_U), width = 0.2, linewidth = 1, position = position_dodge(width = 0.4)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(breaks = c("M", "F"), values = c("#1d3557", "#c1121f")) +
  facet_wrap(~Quantity) +
  theme(
    axis.title = element_blank(),
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "left"
  )

ggsave(
  filename = "alpha.pdf",
  path = IMGFOLDER,
  plot = alpha_plot,
  width = 6, height = 3
)


## NMI ----
NMI_PPE <- array(NA, dim = c(nyears, nyears, 2))
for (t1 in 1:nyears) {
  for (t2 in 1:nyears) {
    NMI_PPE[t1, t2, 1] <- aricode::NMI(PPE_male[, t1], PPE_male[, t2])
    NMI_PPE[t1, t2, 2] <- aricode::NMI(PPE_female[, t1], PPE_female[, t2])
  }
}
dimnames(NMI_PPE) <- list(
  Year1 = seq(2000, 2020, by = 5),
  Year2 = seq(2000, 2020, by = 5),
  Sex = c("NMI - Male", "NMI - Female")
)

NMI_PPE_plot <- NMI_PPE %>%
  reshape2::melt() %>%
  mutate(col_text = if_else(Year1 == Year2, "black", "white")) %>%
  ggplot(aes(x = Year1, y = Year2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_viridis_c("NMI", limits = c(0, 1)) +
  geom_text(aes(label = round(value, 2), color = col_text)) +
  scale_color_manual(values = c("black", "white")) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  facet_wrap(~Sex, ncol = 2) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  ) +
  coord_fixed(ratio = 1)

ggsave(
  filename = "NMI_PPE.pdf",
  path = IMGFOLDER,
  plot = NMI_PPE_plot,
  width = 6, height = 3.5
)


## Plot them together ----
ggsave(
  filename = "PPE_dependence.pdf",
  path = IMGFOLDER,
  plot = gridExtra::arrangeGrob(
    grobs = list(alpha_plot, NMI_PPE_plot), ncol = 2,
    widths = c(1, 1)
  ),
  width = 12, height = 3.5
)

# SELF-HARM ----
# Number of self-harms per year
ghe %>%
  group_by(Year) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup()

# Number of self-harms per year and sex
ghe %>%
  group_by(Year, Sex) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% pivot_wider(id_cols = Year, names_from = Sex, values_from = Selfharm)

# Number of self-harms per year and sex
ghe %>%
  group_by(Year, Sex) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% pivot_wider(id_cols = Year, names_from = Sex, values_from = Selfharm)


## Plot by sex ----
maxselfharm = ghe %>% 
  group_by(Year, Age, Sex) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% pull(Selfharm) %>% max

plt_nmbr_selfharm_sex = 
  ghe %>% 
  group_by(Year, Age, Sex) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm")) %>%
  ungroup() %>% 
  ggplot(aes(x = as.integer(Age), y = Selfharm)) +
  geom_rect(aes(xmin = 3.5, xmax = 5.5, ymin = 0, ymax = maxselfharm,
                fill = "Teen"), alpha = 0.025) +
  geom_rect(aes(xmin = 5.5, xmax = 7.5, ymin = 0, ymax = maxselfharm, 
                fill = "YAd"), alpha = 0.025) +
  geom_rect(aes(xmin = 7.5, xmax = 10.5, ymin = 0, ymax = maxselfharm,
                fill = "Ad"), alpha = 0.025) +
  geom_line(aes(color = Sex), linewidth = 0.85) + 
  geom_point(aes(color = Sex), size = 1) + 
  scale_fill_manual(values = c("grey90", "grey75", "grey55"),
                    name = "Age",
                    breaks = c("Teen", "YAd", "Ad"),
                    labels = c("Adolescence", "Early Adulthood", "Adulthood")) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), 
                     labels = ages[seq(1, nages, by = 3)]) +
  scale_y_continuous(limits = c(0, maxselfharm)) +
  scale_color_manual(values = c("#1d3557", "#c1121f"),breaks = c("M", "F")) +
  guides(color = guide_legend(override.aes = list(linewidth = 3))) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 3),
                             title.position = "left")) +
  facet_wrap(~Year, ncol = 5) + 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.size = unit(0.2, "in"),
        legend.position = "bottom",
        text = element_text(size = 16))

plt_nmbr_selfharm_sex

width = 12; height = 3; zoom = 1
ggsave(filename = "selfharm_nmbr_sex.pdf",
       plot = plt_nmbr_selfharm_sex,
       path = IMGFOLDER,
       width = width*zoom, height = height*zoom)


## Plot by regions ----
plt_nmbr_selfharm_regions = 
  ghe %>% 
  group_by(Year, Age, RegionN) %>%
  summarise(Selfharm = sum(CauseS == "Self-harm"),
            n = n(),
            Selfharm.prop = Selfharm/n) %>%
  ungroup() %>% 
  ggplot(aes(x = as.integer(Age), y = Selfharm.prop)) +
  geom_rect(aes(xmin = 3.5, xmax = 5.5, ymin = 0, ymax = 1,
                fill = "Teen"), alpha = 0.025) +
  geom_rect(aes(xmin = 5.5, xmax = 7.5, ymin = 0, ymax = 1, 
                fill = "YAd"), alpha = 0.025) +
  geom_rect(aes(xmin = 7.5, xmax = 10.5, ymin = 0, ymax = 1,
                fill = "Ad"), alpha = 0.025) +
  geom_line(aes(color = RegionN), linewidth = 0.85) + 
  geom_point(aes(color = RegionN), size = 1) + 
  scale_fill_manual(values = c("grey90", "grey75", "grey55"),
                    name = "Age",
                    breaks = c("Teen", "YAd", "Ad"),
                    labels = c("Adolescence", "Early Adulthood", "Adulthood")) +
  scale_x_continuous(breaks = seq(1, nages, by = 3), 
                     labels = ages[seq(1, nages, by = 3)]) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = palette_regions) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)),
         fill = "none") +
  facet_wrap(~Year, ncol = 5) + 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.size = unit(0.2, "in"),
        legend.position = "bottom",
        text = element_text(size = 16))

plt_nmbr_selfharm_regions

width = 12; height = 3.25; zoom = 1
ggsave(filename = "selfharm_nmbr_regions.pdf",
       plot = plt_nmbr_selfharm_regions,
       path = IMGFOLDER,
       width = width*zoom, height = height*zoom)
